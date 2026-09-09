////////////////////////////////////////////////////////////////////////////////////////////
// ORIGINAL SCRIPT DEVELOPED BY **GIANFRANCO INGRATTA** (York University, ingratta@yorku.ca)
// WITH SUPPORT FROM SUPERVISOR BRUCE HOWARD (York University, blhoward@yorku.ca)
////////////////////////////////////////////////////////////////////////////////////////////

// Script updated by J. L. Barrow (UMN, jbarrow@umn.edu)
// with **much help** from Anthropic Claude Fable 5 High effort

// =====================================================================
//  get_data.cpp
//
//  Fetch per-spill NuMI beam-quality data for a given time window and
//  write it out as (1) a ROOT TTree ("beam") and (2) a tidy CSV.
//
//  Intended use: associating DUNE ND 2x2 demonstrator triggers with the
//  NuMI beam conditions of the corresponding spill (protons on target,
//  horn current and polarity, beam position and width at the target).
//
//  ---------------------------------------------------------------
//  ACNET vs. the IFBeam database — which is which?
//  ---------------------------------------------------------------
//  The device names used below ("E:TRTGTD", "E:NSLINA", ...) are ACNET
//  device identifiers. ACNET is Fermilab's accelerator control network;
//  every instrument in the accelerator complex (toroids, BPMs, power
//  supplies, ...) is addressed by such a name, with the leading letter
//  ("E:") denoting the ACNET device class/area (external-beamline
//  devices, which includes the NuMI primary proton line).
//
//  This program, however, never talks to ACNET directly. It queries the
//  IFBeam database (IFDB, "Intensity Frontier Beam Database"): an
//  archive maintained for the IF experiments (MINOS, MINERvA, NOvA,
//  DUNE/2x2, ...) into which a curated subset of ACNET device readings
//  is copied, timestamped by accelerator timing (TCLK) events. IFDB is
//  exposed through a public, unauthenticated HTTPS REST API — no
//  Kerberos, VPN, or Fermilab credentials are needed. So: the
//  *variables* are ACNET devices; the *source we read them from* is the
//  IFBeam database's archive of those devices.
//
//  Each fetched reading is tagged with the TCLK event on which it was
//  recorded. We request event "$A9" (query string "e=e,a9"), the NuMI
//  extraction timing event, so that we obtain exactly one reading per
//  NuMI spill (~1 per 1.2-2 s) rather than the continuous datalogger
//  stream.
//
//  ---------------------------------------------------------------
//  References (publicly accessible unless noted)
//  ---------------------------------------------------------------
//  [1] P. Adamson et al., "The NuMI Neutrino Beam",
//      Nucl. Instrum. Meth. A 806 (2016) 279-306, arXiv:1507.06690,
//      https://arxiv.org/abs/1507.06690
//      Sec. 4 describes the primary-beam instrumentation quoted below:
//      the two POT toroids (TOR101 upstream in the NuMI line, TORTGT
//      ~10 m before the target), the capacitative beam position
//      monitors (BPMs), the profile-monitor multiwires, and the
//      focusing horns with forward (FHC, neutrino) / reverse (RHC,
//      antineutrino) current operation.
//  [2] NOvA, IFDBSpillInfo_module.cc — the de-facto reference
//      implementation for turning these IFDB devices into per-spill
//      beam-quality quantities. Public mirror:
//      https://github.com/novaexperiment/novasoft/blob/main/IFDBSpillInfo/IFDBSpillInfo_module.cc
//      (historical link: https://cdcvs.fnal.gov/redmine/projects/novaart/repository/entry/trunk/IFDBSpillInfo/IFDBSpillInfo_module.cc#L685
//      — Fermilab's Redmine now requires SSO). The per-device notes and
//      the horn-current calibration used below are taken from it.
//  [3] IFBeam REST API syntax, "DataAccessSyntax" wiki of the
//      "ifbeamdata" project,
//      https://cdcvs.fnal.gov/redmine/projects/ifbeamdata/wiki/DataAccessSyntax
//      (now behind Fermilab SSO; snapshots exist on web.archive.org).
//      Documents the query parameters used in createUrl(): v=device,
//      e=event, t0/t1=time window, f=output format.
//
//  ---------------------------------------------------------------
//  Beam-quality ("goodbeam") cuts: NOvA's criteria vs. this code
//  ---------------------------------------------------------------
//  NOvA's per-spill beam-quality selection is defined in [2] (the
//  goodbeam block at the end of IFDBSpillInfo_module.cc) with the
//  numerical values set in the accompanying IFDBSpillInfo.fcl
//  (standard_ifdbspillinfo). A spill is "goodbeam" only if ALL of the
//  following hold:
//
//    1. POT:          spillpot > 2.00e12 protons (MinPOTCut, quoted in
//                     units of 1e12). POT is E:TRTGTD, falling back to
//                     E:TR101D when TRTGTD is missing or reads < 0.02;
//                     negative toroid readings are clamped to 0.
//    2. Horn current: -202 kA < I_horn < -196.4 kA (Min/MaxHornICut,
//                     FHC sign convention), where I_horn is the
//                     calibrated sum of E:NSLINA-D computed exactly as
//                     in this file's NSLIN block. Spills with
//                     |I_horn| < 1 kA (HornICut0HC) are flagged
//                     horn-off rather than good.
//    3. Position x:   -2.00 mm < posx < +2.00 mm (Min/MaxPosXCut).
//    4. Position y:   -2.00 mm < posy < +2.00 mm (Min/MaxPosYCut;
//                     older NOvA run epochs used (-4.50, 0.00) and
//                     (-3.50, +0.50) — the module selects by run
//                     number). NOTE: NOvA's posx/posy are NOT the raw
//                     E:H/VPTGT readbacks: they are the
//                     intensity-weighted per-batch positions (BPM
//                     array element 0, the auto-tune average, skipped;
//                     batches with zero BPM intensity skipped)
//                     linearly extrapolated from the 121 and TGT BPM
//                     stations to the target z (BpmProjection in [2]).
//    5. Width:        0.57 mm < widthx, widthy < 1.58 mm
//                     (Min/MaxWidthXCut/YCut; 1.88 mm maximum in later
//                     epochs), where the widths are Gaussian sigmas
//                     fitted to the E:MTGTDS[] multiwire profile
//                     (horizontal wires at array indices 103-150,
//                     vertical at 151-198; ProfileProjection in [2]).
//    6. Timing:       |trigger time - IFDB spill time| < 0.5 s
//                     (MaxDeltaTCut = 0.5e9 ns) — the detector trigger
//                     must actually match the spill record.
//
//  What THIS code applies: the "beam" tree and main CSV remain
//  deliberately UNCUT — every spill returned by IFBeam is written out
//  so downstream studies can see off-nominal beam. In addition, a
//  parallel "quality" TTree (and a *_quality.csv), entry-aligned with
//  "beam" so it can be attached as a friend tree, EVALUATES cuts 1-5
//  per spill exactly as NOvA does — the BpmProjection / BpmAtTarget /
//  ProfileProjection / GetGaussFit algorithms of [2] are ported below
//  — and stores the derived quantities (spillpot, hornI, posx, posy,
//  widthx, widthy), the six individual pass flags, and the combined
//  goodbeam verdict. Nothing is filtered: failing spills are written
//  with goodbeam = 0. Cut 6 (trigger timing) is NOT evaluated — it
//  compares a detector trigger time against the spill record, and
//  this standalone fetcher has no detector trigger stream — so
//  "goodbeam" here means cuts 1-5. The only other selections made
//  here are data-hygiene:
//    - readings from different devices are associated to the same
//      spill only if their timestamps agree within beamMatchDT = 0.5 s
//      (the same coincidence half-width as NOvA's timing cut 6);
//    - the derived E:NSLIN horn current is computed only for spills
//      where all four striplines report a nonzero, time-matched
//      reading (a presence requirement, not a quality cut);
//    - unmatched devices appear as NaN (scalars) or empty vectors
//      (arrays) in the TTree — nothing is silently dropped.
//
//  The plotting scripts in plots/ (beam_gif.py, beam_batches_gif.py,
//  beam_slosh.C) additionally apply, following BpmProjection in [2]:
//    - BPM array element 0 (auto-tune average) excluded;
//    - exact-zero BPM readings excluded (the acquisition's bad-batch
//      convention);
//    - horizontal/vertical readings paired within +-0.5 s;
//  and they DRAW the +-2 mm NOvA position box (cuts 3-4) on the
//  figures for illustration only — it is not applied as a filter, and
//  it is drawn against the raw TGT-station readback rather than
//  NOvA's extrapolated-to-target position (see the NOTE under cut 4).
// =====================================================================

#include <curl/curl.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>

#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TF1.h>

using json = nlohmann::json;

// One device's archive for the queried window:
//   key   = spill time, seconds since the Unix epoch (IFBeam "clock"),
//   value = the reading(s) on that spill — a single number for scalar
//           devices (toroids, horn striplines), an array for the
//           BPM ("[]") and multiwire devices.
using BeamInfo = std::map<double, std::vector<double>>;
using DeviceMap = std::map<std::string, BeamInfo>;

// Half-width (seconds) of the coincidence window used to decide that
// readings from two different devices belong to the same spill. Spills
// are >~1.2 s apart, so +-0.5 s matches at most one spill.
double beamMatchDT = 0.5;

// ---------------------------------------------------------------------
// URL construction for the IFBeam REST API [3].
//
// Endpoint: https://dbdata3vm.fnal.gov:9443/ifbeam/data/data
// Query parameters:
//   v  = ACNET device name; a trailing "[]" requests the full array
//        readback rather than only the first element,
//   e  = TCLK event on which the reading was archived; "e,a9" is the
//        NuMI extraction event $A9 (one reading per spill),
//   t0 = window start, ISO-8601 with UTC offset,
//   t1 = window end,
//   f  = response format ("json" here; csv/xml also exist).
// ---------------------------------------------------------------------
std::string createUrl(const std::string& potDevice,
                       const std::string& min_time_iso,
                       const std::string& max_time_iso) {
    // Device names (e.g. "E:HPTGT[]") contain characters that need to be
    // percent-encoded to be valid in a query string.
    CURL* curl = curl_easy_init();
    if (!curl) {
        throw std::runtime_error("createUrl: failed to initialize curl for escaping");
    }
    char* escaped = curl_easy_escape(curl, potDevice.c_str(), static_cast<int>(potDevice.size()));
    std::string encodedDevice = escaped ? escaped : potDevice;
    curl_free(escaped);
    curl_easy_cleanup(curl);

    std::ostringstream url_stream;
    url_stream << "https://dbdata3vm.fnal.gov:9443/ifbeam/data/data?v=" << encodedDevice
               << "&e=" << "e,a9"
               << "&t0=" << min_time_iso
               << "&t1=" << max_time_iso
               << "&f=json";
    return url_stream.str();
}

// ---------------------------------------------------------------------
// Unit-string -> scale factor.
//
// IFBeam stores each reading's engineering unit as a string, and for
// the toroids that string *is* the power of ten, e.g. "E12": a TRTGTD
// value of 49.77 with units "E12" means 4.977e13 protons. This helper
// turns "E<n>" into 10^n and leaves genuine units (e.g. "mm") at 1.0.
// ---------------------------------------------------------------------
double unitToFactor(const std::string& unit) { //because the exponent part in IFBeam is stored as a string under "value" :/
      std::regex re("E(\\d+)");
      std::smatch match;

      if (std::regex_search(unit, match, re) && match.size() > 1) {
          int exponent = std::stoi(match.str(1));
          return std::pow(10, exponent);
      } else {
          return 1.0;
      }
  }

// ---------------------------------------------------------------------
// Epoch (seconds, possibly with a fractional part) -> ISO8601 conversion
// ---------------------------------------------------------------------

std::string toISO8601(double time_sec) {
    time_t seconds = static_cast<time_t>(time_sec);
    long long ms = static_cast<long long>(std::llround((time_sec - static_cast<double>(seconds)) * 1000.0));

    // std::localtime is NOT thread-safe (it uses a shared static buffer).
    // We call it only once and reuse the result to avoid both the
    // thread-safety issue and redundant calls.
    std::tm local_tm{};
#if defined(_WIN32)
    localtime_s(&local_tm, &seconds);
#else
    localtime_r(&seconds, &local_tm);
#endif

    std::ostringstream ss;
    ss << std::put_time(&local_tm, "%Y-%m-%dT%H:%M:%S");
    if (ms > 0) {
        ss << '.' << std::setfill('0') << std::setw(3) << ms;
    }

    long gmtoff = local_tm.tm_gmtoff;
    ss << (gmtoff >= 0 ? '+' : '-')
       << std::setfill('0') << std::setw(2) << std::abs(gmtoff) / 3600
       << ':' << std::setfill('0') << std::setw(2) << (std::abs(gmtoff) % 3600) / 60;

    return ss.str();
}

// ---------------------------------------------------------------------
// curl write callback
// ---------------------------------------------------------------------
static size_t WriteCallback(void* contents, size_t size, size_t nmemb, std::string* userp) {
    userp->append(static_cast<char*>(contents), size * nmemb);
    return size * nmemb;
}

// ---------------------------------------------------------------------
// Query one device from the IFBeam database and unpack the JSON.
//
// The response is {"rows": [...]}, one row per spill, each row holding:
//   "clock" — the reading's timestamp in *milliseconds* since the Unix
//             epoch (converted to seconds here),
//   "units" — engineering-unit string, see unitToFactor(),
//   "c"     — the channel array; each element {"v": <number>} is one
//             reading (1 element for scalar devices; 7 for the BPM
//             arrays; 216 for the MTGTDS multiwire).
// ---------------------------------------------------------------------
BeamInfo retrieveInfoFromDataBase(const std::string& url) {
    constexpr double ms_to_s = 1e-3;
    BeamInfo data;

    CURL* curl = curl_easy_init();
    if (!curl) {
        throw std::runtime_error("retrieveInfoFromDataBase: curl_easy_init() failed");
    }

    std::string readBuffer;
    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &readBuffer);
    curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
    curl_easy_setopt(curl, CURLOPT_TIMEOUT, 30L);

    CURLcode res = curl_easy_perform(curl);
    curl_easy_cleanup(curl);

    if (res != CURLE_OK) {
        // std::abort() used to kill the entire process for a single network
        // error. Better to throw an exception: the caller decides whether
        // to retry, skip the device, or terminate.
        throw std::runtime_error("curl_easy_perform() failed: " +
                                  std::string(curl_easy_strerror(res)));
    }

    try {
        auto json_data = json::parse(readBuffer);
        for (const auto& spill : json_data["rows"]) {
            double time = spill["clock"].get<double>() * ms_to_s;
            std::string unit = spill.value("units", "");
            std::vector<double> values;

            if (spill.contains("c")) {
                for (const auto& c_entry : spill["c"]) {
                    if (c_entry.contains("v") && c_entry["v"].is_number()) {
                        values.push_back(c_entry["v"].get<double>() * unitToFactor(unit));
                    }
                }
            }
            data[time] = std::move(values);
        }
    } catch (const json::exception& e) {
        throw std::runtime_error("JSON parsing error for url '" + url + "': " + e.what());
    }

    return data;
}

// ---------------------------------------------------------------------
// Find the entry of `info` whose spill time is within `dt` seconds of `t`,
// or nullptr if there is none.
// ---------------------------------------------------------------------
const std::vector<double>* findMatch(const BeamInfo& info, double t, double dt) {
    auto it = info.lower_bound(t - dt);
    if (it != info.end() && std::fabs(it->first - t) <= dt) return &it->second;
    return nullptr;
}

// =====================================================================
// NOvA-style beam-quality computation, ported from IFDBSpillInfo [2].
// =====================================================================

// BPM station z-positions [feet], surveyed, origin between the two
// Budal monitors of the NOvA target; only ratios enter the linear
// extrapolation, so the unit is irrelevant. SWIC wire pitch in mm.
constexpr double kZhp121 = -68.04458, kZvp121 = -66.99283,
                 kZhptgt = -31.25508, kZvptgt = -30.16533, kZtarg = 0.0;
constexpr double kWireSpacing = 0.5;   // mm

// Port of extrapolate_position [2]: linear extrapolation of a
// transverse position measured at z1 and z2 to z3.
double extrapolatePosition(double t1, double z1, double t2, double z2, double z3) {
    return t1 + (t2 - t1) * (z3 - z1) / (z2 - z1);
}

// Port of BpmProjection + BpmAtTarget [2]: intensity-weighted mean of
// the per-batch BPM positions extrapolated to the target z. Array
// element 0 (the auto-tune average) is skipped, as are batches whose
// BPM intensity reads zero. Positions in mm. Returns false when no
// valid batch remains. NOTE: NOvA's BpmAtTarget reads the *horizontal*
// intensity into both weights (IY = BpmIntX[ind], an apparent slip);
// here each plane is weighted by its own intensity — the two readbacks
// are nearly identical, so the difference is negligible.
bool computeBpmPosition(const std::vector<double>& vp121, const std::vector<double>& hp121,
                        const std::vector<double>& vptgt, const std::vector<double>& hptgt,
                        const std::vector<double>& vitgt, const std::vector<double>& hitgt,
                        double& posx, double& posy) {
    size_t n = std::min({vp121.size(), hp121.size(), vptgt.size(),
                         hptgt.size(), vitgt.size(), hitgt.size()});
    if (n == 0) return false;
    size_t start = (n > 1) ? 1 : 0;            // skip auto-tune element 0
    double IXtot = 0, IYtot = 0, ix = 0, iy = 0;
    for (size_t ind = start; ind < n; ++ind) {
        double IX = hitgt[ind], IY = vitgt[ind];
        if (IX == 0.0 || IY == 0.0) continue;  // bad batch reading
        if (IX < 0) IX = 0;
        if (IY < 0) IY = 0;
        double xb = extrapolatePosition(hp121[ind], kZhp121, hptgt[ind], kZhptgt, kZtarg);
        double yb = extrapolatePosition(vp121[ind], kZvp121, vptgt[ind], kZvptgt, kZtarg);
        IXtot += IX; IYtot += IY;
        ix += IX * xb; iy += IY * yb;
    }
    if (IXtot <= 0.0 || IYtot <= 0.0) return false;
    posx = ix / IXtot;
    posy = iy / IYtot;
    return true;
}

// Port of GetStats [2]: charge statistics of the inverted wire profile,
// used as fit start values. As in NOvA, the graph is shifted up by its
// most negative value as a side effect (so the subsequent fit runs on
// the shifted profile, and the offset start value falls back to -3
// because no negative readings remain).
double profileStats(TGraph& prof, double& mean, double& stdev) {
    mean = 0.0; stdev = 0.0;
    double minprof = 1e6;
    for (int i = 0; i < prof.GetN(); ++i) minprof = std::min(minprof, prof.GetY()[i]);
    double qtot = 0, qx = 0;
    for (int i = 0; i < prof.GetN(); ++i) {
        prof.SetPoint(i, prof.GetX()[i], prof.GetY()[i] - minprof);
        qtot += prof.GetY()[i];
        qx += prof.GetY()[i] * prof.GetX()[i];
    }
    if (qtot == 0.0) return qtot;
    mean = qx / qtot;
    double var = 0;
    for (int i = 0; i < prof.GetN(); ++i) {
        double d = prof.GetX()[i] - mean;
        var += prof.GetY()[i] * d * d;
    }
    var /= qtot;
    if (var > 0.0) stdev = std::sqrt(var);
    return qtot;
}

// Port of ProfileProjection + GetGaussFit [2]: Gaussian fit to one
// plane's 48 multiwire channels (raw voltages [V]; NOvA works in mV,
// hence the x1000). Wires sit at (ch - 24.5) * 0.5 mm; the profile is
// inverted so the peak is positive; dead (exactly zero) channels are
// dropped. Returns the fitted sigma [mm], NaN when NOvA's
// preconditions or physicality checks fail.
double profileGaussWidth(const std::vector<double>& wires) {
    const double bad = std::numeric_limits<double>::quiet_NaN();
    if (wires.size() < 48) return bad;
    TGraph prof;
    for (int ch = 1; ch <= 48; ++ch) {
        double xpos = (ch - 24.5) * kWireSpacing;   // mm
        double qmv = wires[ch - 1] * 1000.0;        // V -> mV
        prof.SetPoint(prof.GetN(), xpos, -qmv);     // invert
    }
    for (int i = 0; i < prof.GetN(); ++i)
        if (prof.GetY()[i] == 0.0) { prof.RemovePoint(i); --i; }

    double statMean = 0, statStdev = 0;
    double area = profileStats(prof, statMean, statStdev);   // shifts profile

    // NOvA's "reasonable fit" preconditions: enough live wires, enough
    // charge (mV), stats mean within 10 mm, stats stdev within 6 mm.
    if (!(prof.GetN() > 5 && std::abs(area) > 100.0 &&
          std::abs(statMean) < 20.0 * kWireSpacing &&
          std::abs(statStdev) < 12.0 * kWireSpacing))
        return bad;

    TF1 gausF("gausF",
        "[3]+([0]*exp((-1*(x-[1])*(x-[1]))/(2*[2]*[2]))/(sqrt(2*3.142)*[2]))");
    gausF.SetParameters(area, statMean, statStdev / 2.0, -3.0);
    if (prof.Fit(&gausF, "q", "", -24 * kWireSpacing, 24 * kWireSpacing) != 0)
        return bad;

    // NOvA's physicality checks on the fit result
    if (gausF.GetParameter(0) < 10.0) return bad;
    if (std::abs(gausF.GetParameter(1)) > 24 * kWireSpacing) return bad;
    if (gausF.GetParameter(2) > 24 * kWireSpacing || gausF.GetParameter(2) < 0.0) return bad;
    return gausF.GetParameter(2);   // sigma [mm]
}

// ---------------------------------------------------------------------
// Main entry point (ROOT macro convention: function name == file name).
//
//   root -l -b -q 'get_data.cpp("2024-07-12T00:01:00-05:00","2024-07-12T01:01:00-05:00")'
//
// Times are ISO-8601 with an explicit UTC offset (-05:00 = CDT).
// ---------------------------------------------------------------------
int get_data(const std::string& min_time_iso = "2024-07-12T00:01:00-05:00",
             const std::string& max_time_iso = "2024-07-12T01:01:00-05:00") {

    // -----------------------------------------------------------------
    // The ACNET devices to pull from the IFBeam archive. All are real
    // ACNET device names except E:NSLIN, which is a pseudo-device
    // computed below. Per-device notes follow [1] Sec. 4 and the inline
    // documentation of [2].
    // -----------------------------------------------------------------
    DeviceMap deviceMap = {
        //_________POT (protons on target, per spill)_________
        // Two toroidal beam-current transformers encircle the primary
        // proton beam [1]: TOR101 sits in the upstream NuMI beamline,
        // TORTGT ~10 m before the target. The trailing "D" selects the
        // $A9-datalogged readback of each toroid; [2] describes
        // E:TRTGTD as "intensity (POT) on NuMI target toroid" and
        // E:TR101D as "intensity (POT) on NuMI toroid 101".
        // IFBeam units are "E12", i.e. readings come back in units of
        // 1e12 protons (unitToFactor() applies the factor, so values
        // here are absolute protons/spill, ~5e13 at nominal intensity).
        // [2] uses TRTGTD as the primary POT device and falls back to
        // TR101D when TRTGTD is missing or unphysically small.
        {"E:TRTGTD", {}},   // POT at target toroid TORTGT (primary POT device)
        {"E:TOR101", {}},   // upstream toroid TOR101, direct readback
        {"E:TR101D", {}},   // upstream toroid TOR101, datalogged readback (POT backup)

        //_________Horn current_______________________________
        // The NuMI focusing-horn current (~200 kA nominal [1]) is
        // delivered through parallel stripline conductors; NSLINA-D are
        // the four stripline current transductor readbacks (each ~1/4
        // of the total, in kA). [2] documents them as "current on NuMI
        // strip line A/B/C/D"; the calibrated sum (see below) gives the
        // total horn current.
        {"E:NSLINA", {}},   // stripline A current [kA]
        {"E:NSLINB", {}},   // stripline B current [kA]
        {"E:NSLINC", {}},   // stripline C current [kA]
        {"E:NSLIND", {}},   // stripline D current [kA]
        // Pseudo-device: NOT fetched from IFBeam. Filled below with the
        // calibrated total horn current computed from NSLINA-D.
        {"E:NSLIN", {}},    // total horn current [kA] (computed here)

        //________Horn polarity_______________________________
        // Horn direction readback: distinguishes forward horn current
        // (FHC, focuses pi+ -> neutrino beam) from reverse horn current
        // (RHC -> antineutrino beam) [1]. NOTE: no public document
        // describing this device's value convention was found ([2] does
        // not use it; NOvA takes FHC/RHC from its run database
        // instead), so validate the readback against a known beam
        // configuration before relying on it.
        {"E:HRNDIR", {}},

        //__________Beam position (BPMs)______________________
        // Capacitative beam position monitors in the primary line [1]:
        // station "121" is upstream in the NuMI line, "TGT" is the last
        // station before the target. H/V prefix = horizontal/vertical
        // plane; "P" = position [mm], "I" = the intensity seen by the
        // BPM (used to weight/validate the position measurements).
        // The "[]" suffix requests the full array: 7 elements per
        // spill. Per [2]: element 0 is "an average value of a few
        // batches (often 2,3,4) used for auto-tuning, so we should
        // disregard it"; elements 1-6 are the per-Booster-batch
        // readings (a NuMI spill is 6 Booster batches). [2] linearly
        // extrapolates the 121- and TGT-station positions batch-by-
        // batch to the target z to obtain the beam spot position.
        {"E:HPTGT[]", {}},  // horizontal position at target BPM [mm], 7 elements
        {"E:HITGT[]", {}},  // intensity, horizontal plane, target BPM
        {"E:HP121[]", {}},  // horizontal position at station-121 BPM [mm]
        {"E:VPTGT[]", {}},  // vertical position at target BPM [mm]
        {"E:VITGT[]", {}},  // intensity, vertical plane, target BPM
        {"E:VP121[]", {}},  // vertical position at station-121 BPM [mm]

        //__________Beam width (profile monitor)______________
        // Multiwire profile monitor at the target ("MTGTDS") [1]. The
        // 216-element array readback contains, per [2]
        // (ProfileProjection): the 48 horizontal-wire charges at
        // indices [103..150] and the 48 vertical-wire charges at
        // [151..198] (remaining elements are header/housekeeping
        // words). A Gaussian fit to each plane's wire profile yields
        // the beam centroid and width (sigma) in mm at the target.
        {"E:MTGTDS[]", {}}
    };

    // -----------------------------------------------------------------
    // Fetch every device's archive for the requested window.
    // Failures (network, device missing from the archive) skip that
    // device and continue, mirroring the tolerant behavior of [2].
    // -----------------------------------------------------------------
    for (auto& [deviceName, beamInfo] : deviceMap) {
        beamInfo.clear();
        std::string url_device = createUrl(deviceName, min_time_iso, max_time_iso);
        std::cout << url_device << "\n";
        try {
            beamInfo = retrieveInfoFromDataBase(url_device);
        } catch (const std::exception& e) {
            std::cerr << "Error on device '" << deviceName << "': " << e.what() << "\n";
            // continue with the other devices instead of crashing everything
        }
        std::cout << deviceName << ": " << beamInfo.size() << " spills\n";
    }

    // -----------------------------------------------------------------
    // Total horn current from the four stripline readbacks.
    //
    // Calibration (per-stripline offset and gain) taken verbatim from
    // NOvA's IFDBSpillInfo [2], where it is attributed to
    // "Updated numbers from Jim Hylen, email Oct 4, 2013":
    //   I_horn = (E:NSLINA-(+0.01))/0.9951 + (E:NSLINB-(-0.14))/0.9957
    //          + (E:NSLINC-(-0.05))/0.9965 + (E:NSLIND-(-0.07))/0.9945
    // Striplines B-D are matched to each NSLINA spill within
    // beamMatchDT; spills where any stripline reads exactly 0 are
    // skipped. Result is stored under the pseudo-device E:NSLIN, in kA
    // (~ -200 kA at nominal FHC running in the 2024 data).
    // -----------------------------------------------------------------
    double currentA=0.0, currentB=0.0, currentC=0.0, currentD=0.0;

    for(const auto& pairA : deviceMap["E:NSLINA"]) {
      auto timeA = pairA.first;
      currentA = pairA.second.at(0);

      if(currentA == 0.) continue;

      for(const auto& pairB : deviceMap["E:NSLINB"]) {
        auto timeB = pairB.first;
        if (abs(timeA - timeB) <= beamMatchDT) {
          currentB = pairB.second.at(0);
          break;
        }
      }

      if(currentB == 0.) continue;

      for(const auto& pairC : deviceMap["E:NSLINC"]) {
        auto timeC = pairC.first;
        if (abs(timeA - timeC) <= beamMatchDT)
        {
          currentC = pairC.second.at(0);
          break;
        }
      }

      if(currentC == 0.) continue;

      for(const auto& pairD : deviceMap["E:NSLIND"]) {
        auto timeD = pairD.first;
        if (abs(timeA - timeD) <= beamMatchDT){
          currentD = pairD.second.at(0);
          break;
        }
      }

      if(currentD == 0.) continue;

      deviceMap["E:NSLIN"][timeA] = {(currentA-0.01)/0.9951 + (currentB+0.14)/0.9957 + (currentC+0.05)/0.9965 + (currentD+0.07)/0.9945};
    }

    // -----------------------------------------------------------------
    // Output 1: ROOT TTree, one entry per spill.
    // The device with the most spills defines the reference timeline;
    // every other device is matched to it within beamMatchDT. Unmatched
    // scalars are NaN, unmatched arrays are empty vectors.
    // -----------------------------------------------------------------
    std::string refDevice;
    size_t refCount = 0;
    for (const auto& [name, info] : deviceMap)
        if (info.size() > refCount) { refCount = info.size(); refDevice = name; }
    if (refCount == 0) {
        std::cerr << "No data retrieved for any device; no output written.\n";
        return 1;
    }
    std::cout << "Reference timeline: " << refDevice << " (" << refCount << " spills)\n";

    // Filenames carry the queried time range. ':' is illegal in filenames
    // on Windows/OneDrive, so it is stripped from the ISO timestamps:
    // "2024-07-12T00:01:00-05:00" -> "2024-07-12T000100-0500"
    auto sanitizeTime = [](std::string s) {
        s.erase(std::remove(s.begin(), s.end(), ':'), s.end());
        return s;
    };
    const std::string fileStem = "get_data_" + sanitizeTime(min_time_iso) +
                                 "_to_" + sanitizeTime(max_time_iso);
    const std::string rootFileName = fileStem + ".root";
    const std::string csvFileName  = fileStem + ".csv";

    // "E:HPTGT[]" -> branch/column name "HPTGT"
    auto branchName = [](std::string name) {
        if (name.rfind("E:", 0) == 0) name = name.substr(2);
        while (!name.empty() && (name.back() == '[' || name.back() == ']')) name.pop_back();
        return name;
    };

    // Devices whose spills only ever carry one value get a plain double
    // branch; the rest get a std::vector<double> branch.
    std::map<std::string, bool> isScalar;
    for (const auto& [name, info] : deviceMap) {
        size_t maxSize = 0;
        for (const auto& [t, v] : info) maxSize = std::max(maxSize, v.size());
        isScalar[name] = (maxSize <= 1);
    }

    TFile fout(rootFileName.c_str(), "RECREATE");
    TTree tree("beam", "NuMI IFBeam per-spill data");
    double spillTime = 0.0;
    std::string spillISO;
    tree.Branch("time", &spillTime);
    tree.Branch("time_iso", &spillISO);
    std::map<std::string, double> scalarBuf;
    std::map<std::string, std::vector<double>> vectorBuf;
    for (const auto& [name, info] : deviceMap) {
        const std::string bname = branchName(name);
        if (isScalar[name]) tree.Branch(bname.c_str(), &scalarBuf[name]);
        else                tree.Branch(bname.c_str(), &vectorBuf[name]);
    }

    // -----------------------------------------------------------------
    // Parallel NOvA-style "quality" tree, entry-aligned with "beam" so
    // it can be attached as a friend tree. Evaluates goodbeam cuts 1-5
    // (see the "Beam-quality cuts" header comment; cut 6 needs a
    // detector trigger stream and is not evaluated). Nothing is
    // filtered — failing spills are written with goodbeam = 0.
    // -----------------------------------------------------------------
    // NOvA standard_ifdbspillinfo cut values (IFDBSpillInfo.fcl):
    const double kMinPOT = 2.0e12;                        // protons/spill
    const double kMinHornI = -202.0, kMaxHornI = -196.4;  // kA (FHC)
    const double kHornI0HC = 1.0;                         // kA, horn-off flag
    const double kMinPosX = -2.0, kMaxPosX = 2.0;         // mm
    const double kMinPosY = -2.0, kMaxPosY = 2.0;         // mm (standard epoch)
    const double kMinWidth = 0.57, kMaxWidth = 1.58;      // mm (standard epoch)

    TTree qtree("quality", "NOvA-style per-spill goodbeam evaluation (cuts 1-5)");
    const double dNaN = std::numeric_limits<double>::quiet_NaN();
    double q_pot, q_horn, q_posx, q_posy, q_wx, q_wy;
    Bool_t p_pot, p_horn, p_posx, p_posy, p_wx, p_wy, q_good, q_is0hc;
    qtree.Branch("time", &spillTime);
    qtree.Branch("time_iso", &spillISO);
    qtree.Branch("spillpot", &q_pot);      // NOvA POT logic (TRTGTD w/ TR101D fallback)
    qtree.Branch("hornI", &q_horn);        // calibrated NSLIN sum [kA]
    qtree.Branch("posx", &q_posx);         // intensity-weighted, extrapolated to target [mm]
    qtree.Branch("posy", &q_posy);
    qtree.Branch("widthx", &q_wx);         // multiwire Gaussian sigma [mm]
    qtree.Branch("widthy", &q_wy);
    qtree.Branch("is0HC", &q_is0hc);       // |hornI| < 1 kA (horn off)
    qtree.Branch("pass_pot", &p_pot);
    qtree.Branch("pass_horn", &p_horn);
    qtree.Branch("pass_posx", &p_posx);
    qtree.Branch("pass_posy", &p_posy);
    qtree.Branch("pass_widthx", &p_wx);
    qtree.Branch("pass_widthy", &p_wy);
    qtree.Branch("goodbeam", &q_good);

    std::ofstream qcsv(fileStem + "_quality.csv");
    qcsv << "time_s,time_iso,spillpot,hornI,posx,posy,widthx,widthy,"
            "pass_pot,pass_horn,pass_posx,pass_posy,pass_widthx,pass_widthy,goodbeam\n";
    qcsv << std::setprecision(15);
    long long nGood = 0;

    for (const auto& [t, refValues] : deviceMap[refDevice]) {
        spillTime = t;
        spillISO = toISO8601(t);
        for (const auto& [name, info] : deviceMap) {
            const std::vector<double>* match = findMatch(info, t, beamMatchDT);
            if (isScalar[name]) {
                scalarBuf[name] = (match && !match->empty())
                                      ? match->front()
                                      : std::numeric_limits<double>::quiet_NaN();
            } else {
                vectorBuf[name] = match ? *match : std::vector<double>{};
            }
        }
        tree.Fill();

        // ---- NOvA-style quality evaluation for this spill ----
        // 1. POT: TRTGTD, falling back to TR101D when TRTGTD is missing
        //    or reads below 0.02e12; negative readings clamped to 0.
        double trtgtd = scalarBuf["E:TRTGTD"];
        double tr101d = scalarBuf["E:TR101D"];
        q_pot = dNaN;
        if (!std::isnan(trtgtd)) {
            q_pot = trtgtd;
            if (!std::isnan(tr101d) && tr101d != 0.0 && q_pot < 0.02e12) q_pot = tr101d;
        } else if (!std::isnan(tr101d)) {
            q_pot = tr101d;
        }
        if (q_pot < 0.0) q_pot = 0.0;
        p_pot = !std::isnan(q_pot) && q_pot > kMinPOT;

        // 2. Horn current: the calibrated NSLIN sum computed above.
        q_horn = scalarBuf["E:NSLIN"];
        p_horn = q_horn > kMinHornI && q_horn < kMaxHornI;
        q_is0hc = std::fabs(q_horn) < kHornI0HC;

        // 3-4. Position at target: intensity-weighted BPM extrapolation.
        q_posx = dNaN; q_posy = dNaN;
        computeBpmPosition(vectorBuf["E:VP121[]"], vectorBuf["E:HP121[]"],
                           vectorBuf["E:VPTGT[]"], vectorBuf["E:HPTGT[]"],
                           vectorBuf["E:VITGT[]"], vectorBuf["E:HITGT[]"],
                           q_posx, q_posy);
        p_posx = q_posx > kMinPosX && q_posx < kMaxPosX;
        p_posy = q_posy > kMinPosY && q_posy < kMaxPosY;

        // 5. Width: Gaussian sigma of each multiwire plane
        //    (MTGTDS elements 103-150 horizontal, 151-198 vertical).
        const std::vector<double>& mw = vectorBuf["E:MTGTDS[]"];
        q_wx = dNaN; q_wy = dNaN;
        if (mw.size() >= 199) {
            q_wx = profileGaussWidth({mw.begin() + 103, mw.begin() + 151});
            q_wy = profileGaussWidth({mw.begin() + 151, mw.begin() + 199});
        }
        p_wx = q_wx > kMinWidth && q_wx < kMaxWidth;
        p_wy = q_wy > kMinWidth && q_wy < kMaxWidth;

        q_good = p_pot && p_horn && p_posx && p_posy && p_wx && p_wy;
        if (q_good) ++nGood;
        qtree.Fill();
        qcsv << spillTime << ',' << spillISO << ',' << q_pot << ',' << q_horn << ','
             << q_posx << ',' << q_posy << ',' << q_wx << ',' << q_wy << ','
             << int(p_pot) << ',' << int(p_horn) << ',' << int(p_posx) << ','
             << int(p_posy) << ',' << int(p_wx) << ',' << int(p_wy) << ','
             << int(q_good) << '\n';
    }
    tree.Write();
    qtree.Write();
    qcsv.close();
    std::cout << "Wrote " << tree.GetEntries() << " spills to " << rootFileName
              << " (TTrees 'beam' + 'quality')\n";
    std::cout << "goodbeam (cuts 1-5): " << nGood << " / " << qtree.GetEntries()
              << " spills pass; quality CSV: " << fileStem << "_quality.csv\n";
    fout.Close();

    // -----------------------------------------------------------------
    // Output 2: CSV in long/tidy format, one row per (device, spill,
    // array index) — no time matching, exactly what came back per device.
    // -----------------------------------------------------------------
    std::ofstream csv(csvFileName);
    csv << "device,time_s,time_iso,index,value\n";
    csv << std::setprecision(15);
    long long nRows = 0;
    for (const auto& [name, info] : deviceMap) {
        for (const auto& [t, values] : info) {
            const std::string iso = toISO8601(t);
            for (size_t i = 0; i < values.size(); ++i, ++nRows)
                csv << name << ',' << t << ',' << iso << ',' << i << ',' << values[i] << '\n';
        }
    }
    csv.close();
    std::cout << "Wrote " << nRows << " rows to " << csvFileName << "\n";

    return 0;
}

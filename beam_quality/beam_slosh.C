// =====================================================================
//  beam_slosh.C — animated COLZ heatmap of the NuMI beam position.
//
//  Reads the "beam" TTree written by get_data.cpp and renders an
//  animated GIF: each frame is a TH2D filled with the individual
//  batch-by-batch beam positions (BPM array elements 1-6; element 0 is
//  the auto-tune average and zeros are bad readings, both excluded) of
//  10 consecutive spills, drawn with COLZ. Each batch's 121- and
//  TGT-station readings are linearly extrapolated to the target z
//  (NOvA BpmProjection geometry), so the plotted coordinate is the one
//  the goodbeam position cut is defined in. The blob's frame-to-frame
//  motion shows the beam position "sloshing". The red box is the NOvA
//  goodbeam position window (cuts 3-4) that get_data.cpp recorded in
//  the same ROOT file as cut_minPosX/cut_maxPosX/cut_minPosY/
//  cut_maxPosY (TParameter<double>) and cut_epoch (TNamed); for all
//  2x2 running that is NOvA's periods-13+ window, x in (-2, +2) mm and
//  y in (-3.5, +0.5) mm.
//
//  The full catalogue of NOvA's goodbeam criteria, and exactly which
//  of them this pipeline applies, illustrates, or omits, is documented
//  in the "Beam-quality cuts" comment block of get_data.cpp.
//
//  Usage:
//    root -l -b -q beam_slosh.C
// =====================================================================
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TH2D.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TROOT.h>
#include <TSystem.h>
#include <algorithm>
#include <cstdio>
#include <ctime>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

void beam_slosh(const std::string& inFile =
                    "get_data_2024-07-12T000100-0500_to_2024-07-12T010100-0500.root",
                const std::string& outGif =
                    "beam_position_slosh_root_2024-07-12T000100-0500_to_2024-07-12T010100-0500.gif",
                int spillsPerFrame = 10)
{
    // ---------------- read all spills into memory ----------------
    TFile f(inFile.c_str());
    if (f.IsZombie()) { std::cerr << "cannot open " << inFile << "\n"; return; }
    TTreeReader reader("beam", &f);
    TTreeReaderValue<double> time(reader, "time");
    TTreeReaderValue<std::vector<double>> hp(reader, "HPTGT");
    TTreeReaderValue<std::vector<double>> vp(reader, "VPTGT");
    TTreeReaderValue<std::vector<double>> hu(reader, "HP121");
    TTreeReaderValue<std::vector<double>> vu(reader, "VP121");

    // Surveyed BPM station z-positions [feet] and the NOvA linear
    // extrapolation of each batch to the target (z = 0): the same
    // constants and formula as bpmProjection() in get_data.cpp
    // (ported from NOvA's IFDBSpillInfo extrapolate_position /
    // BpmProjection), so the plotted coordinate matches the one the
    // goodbeam position cut acts on.
    constexpr double zHp121 = -68.04458, zVp121 = -66.99283,
                     zHptgt = -31.25508, zVptgt = -30.16533;
    constexpr double fx = (0.0 - zHp121) / (zHptgt - zHp121);
    constexpr double fy = (0.0 - zVp121) / (zVptgt - zVp121);

    struct Spill { double t; std::vector<std::pair<double,double>> pts; };
    std::vector<Spill> spills;
    while (reader.Next()) {
        size_t n = std::min({hp->size(), vp->size(), hu->size(), vu->size()});
        Spill s{*time, {}};
        for (size_t j = 1; j < n; ++j)              // skip auto-tune element 0
            if ((*hp)[j] != 0.0 && (*vp)[j] != 0.0 &&   // zeros = bad readings
                (*hu)[j] != 0.0 && (*vu)[j] != 0.0)
                s.pts.emplace_back((*hu)[j] + ((*hp)[j] - (*hu)[j]) * fx,
                                   (*vu)[j] + ((*vp)[j] - (*vu)[j]) * fy);
        if (!s.pts.empty()) spills.push_back(std::move(s));
    }
    int nFrames = static_cast<int>(spills.size()) / spillsPerFrame;
    std::cout << spills.size() << " spills with BPM data -> " << nFrames << " frames\n";
    if (nFrames == 0) return;

    // ---------------- fixed axis ranges over the whole window ----------------
    double xlo = 1e9, xhi = -1e9, ylo = 1e9, yhi = -1e9;
    for (const auto& s : spills)
        for (const auto& p : s.pts) {
            xlo = std::min(xlo, p.first);  xhi = std::max(xhi, p.first);
            ylo = std::min(ylo, p.second); yhi = std::max(yhi, p.second);
        }
    // The NOvA "goodbeam" position window (cuts 3-4) exactly as
    // get_data.cpp applied it: read back from the cut_* objects it wrote
    // into this file. Files written before that record existed fall back
    // to NOvA's current values (periods 13+, all 2x2 running), with a
    // warning — the window is asymmetric in y since NOvA period 11.
    auto readCut = [&f](const char* name, double fallback) {
        auto* p = dynamic_cast<TParameter<double>*>(f.Get(name));
        if (!p) {
            std::cerr << "warning: " << name << " not in " << f.GetName()
                      << "; using NOvA periods-13+ value " << fallback << "\n";
            return fallback;
        }
        return p->GetVal();
    };
    const double cutX0 = readCut("cut_minPosX", -2.00), cutX1 = readCut("cut_maxPosX", +2.00);
    const double cutY0 = readCut("cut_minPosY", -3.50), cutY1 = readCut("cut_maxPosY", +0.50);
    auto* epochObj = dynamic_cast<TNamed*>(f.Get("cut_epoch"));
    const std::string epoch = epochObj ? epochObj->GetTitle() : "P13+ (assumed)";

    // Never clip: axes span the union of the data and the goodbeam
    // window, so off-target excursions and the cut boundary are always
    // in view.
    xlo = std::min(xlo, cutX0); xhi = std::max(xhi, cutX1);
    ylo = std::min(ylo, cutY0); yhi = std::max(yhi, cutY1);
    double xpad = 0.08 * (xhi - xlo), ypad = 0.08 * (yhi - ylo);
    xlo -= xpad; xhi += xpad; ylo -= ypad; yhi += ypad;

    TH2D h("h", ";Horizontal position at target [mm];Vertical position at target [mm]",
           96, xlo, xhi, 96, ylo, yhi);

    // Pass 1: global maximum occupancy, so the color scale is fixed
    // across frames (a per-frame autoscale would hide intensity changes).
    double gmax = 0;
    for (int fr = 0; fr < nFrames; ++fr) {
        h.Reset();
        for (int k = 0; k < spillsPerFrame; ++k)
            for (const auto& p : spills[fr * spillsPerFrame + k].pts)
                h.Fill(p.first, p.second);
        gmax = std::max(gmax, h.GetMaximum());
    }

    // ---------------- style ----------------
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kViridis);          // sequential, colorblind-safe
    gStyle->SetNumberContours(64);
    gROOT->SetBatch(kTRUE);

    // TCanvas subtracts window decorations from the requested size in
    // batch; pad the request so the image lands on 3000x2000.
    TCanvas c("c", "c", 3004, 2028);
    c.SetCanvasSize(3000, 2000);
    c.SetRightMargin(0.12);
    c.SetLeftMargin(0.09);
    c.SetTopMargin(0.10);
    c.SetBottomMargin(0.09);

    TLatex title; title.SetNDC(); title.SetTextFont(42);
    TLatex stamp; stamp.SetNDC(); stamp.SetTextFont(102);

    // NOvA goodbeam position window (IFDBSpillInfo: goodbeam requires
    // posx and posy inside it, alongside POT/horn/width/timing cuts).
    TBox cutBox(cutX0, cutY0, cutX1, cutY1);
    cutBox.SetFillStyle(0);
    cutBox.SetLineColor(kRed + 1);
    cutBox.SetLineWidth(4);
    cutBox.SetLineStyle(7);
    TLatex cutLab; cutLab.SetTextFont(42); cutLab.SetTextColor(kRed + 1);

    gSystem->Unlink(outGif.c_str());       // Print("+") appends; start clean

    // Pass 2: draw and append each frame ("+8" = 80 ms per frame).
    for (int fr = 0; fr < nFrames; ++fr) {
        h.Reset();
        double t0 = spills[fr * spillsPerFrame].t;
        for (int k = 0; k < spillsPerFrame; ++k)
            for (const auto& p : spills[fr * spillsPerFrame + k].pts)
                h.Fill(p.first, p.second);
        h.SetMinimum(0);
        h.SetMaximum(gmax);
        h.Draw("COLZ");
        cutBox.Draw("l same");
        cutLab.SetTextSize(0.022);
        cutLab.DrawLatex(cutX0, cutY1 + 0.015 * (yhi - ylo),
                         Form("NuMI good-beam position window (NOvA %s): x %+.2f..%+.2f, y %+.2f..%+.2f mm",
                              epoch.c_str(), cutX0, cutX1, cutY0, cutY1));

        title.SetTextSize(0.033);
        title.DrawLatex(0.09, 0.955,
            "NuMI beam position at target  (per-batch, extrapolated to target z, 10 spills/frame)");
        // spill epoch (UTC) -> CDT = UTC-5
        std::time_t tt = static_cast<std::time_t>(t0) - 5 * 3600;
        char buf[64]; std::strftime(buf, sizeof buf, "%Y-%m-%d %H:%M:%S CDT", std::gmtime(&tt));
        stamp.SetTextSize(0.026);
        stamp.DrawLatex(0.09, 0.915,
            Form("%s   spills %d-%d", buf, fr * spillsPerFrame,
                 fr * spillsPerFrame + spillsPerFrame - 1));

        c.Print((outGif + "+8").c_str());
        if (fr % 25 == 0) std::cout << "frame " << fr << "/" << nFrames << "\n" << std::flush;
    }
    std::cout << "done: " << outGif << "\n";
}

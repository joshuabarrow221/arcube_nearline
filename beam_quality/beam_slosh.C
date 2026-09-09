// =====================================================================
//  beam_slosh.C — animated COLZ heatmap of the NuMI beam position.
//
//  Reads the "beam" TTree written by get_data.cpp and renders an
//  animated GIF: each frame is a TH2D filled with the individual
//  batch-by-batch target-BPM readings (E:HPTGT/E:VPTGT array elements
//  1-6; element 0 is the auto-tune average and zeros are bad readings,
//  both excluded) of 10 consecutive spills, drawn with COLZ. The blob's
//  frame-to-frame motion shows the beam position "sloshing".
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
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
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

    struct Spill { double t; std::vector<std::pair<double,double>> pts; };
    std::vector<Spill> spills;
    while (reader.Next()) {
        size_t n = std::min(hp->size(), vp->size());
        Spill s{*time, {}};
        for (size_t j = 1; j < n; ++j)              // skip auto-tune element 0
            if ((*hp)[j] != 0.0 && (*vp)[j] != 0.0) // zeros = bad batch readings
                s.pts.emplace_back((*hp)[j], (*vp)[j]);
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
    // Never clip: axes span the union of the data and the NOvA
    // "goodbeam" position box (IFDBSpillInfo.fcl: Min/MaxPosXCut,
    // Min/MaxPosYCut = +-2 mm at target), so off-target excursions and
    // the cut boundary are always in view.
    const double cutLo = -2.0, cutHi = +2.0;   // mm
    xlo = std::min(xlo, cutLo); xhi = std::max(xhi, cutHi);
    ylo = std::min(ylo, cutLo); yhi = std::max(yhi, cutHi);
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

    // NOvA goodbeam position box cut (IFDBSpillInfo: goodbeam requires
    // posx and posy inside +-2 mm, alongside POT/horn/width/timing cuts).
    TBox cutBox(cutLo, cutLo, cutHi, cutHi);
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
        cutLab.DrawLatex(cutLo, cutHi + 0.015 * (yhi - ylo),
                         "NOvA goodbeam position cut (#pm2 mm)");

        title.SetTextSize(0.033);
        title.DrawLatex(0.09, 0.955,
            "NuMI beam position at target  (E:HPTGT / E:VPTGT, batch-by-batch, 10 spills per frame)");
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

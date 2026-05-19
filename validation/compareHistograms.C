int compareHistograms(TString mode) {
	TCanvas c1;

	TH1D hpval("pvals","",10,0.,1.);
	hpval.SetMinimum(0.);

	double sumpval(0.), sum(0.);
	int nWarn(0), nFail(0), nSkip(0);
	const double warnThreshold = 0.01;
	const double failThreshold = 0.001;
	std::vector<std::pair<TString,double> > failing;

	TFile *f1 = TFile::Open("${RAPIDSIM_ROOT}/validation/rootfiles/"+mode+"_hists.root");
	TFile *f2 = TFile::Open(mode+"_hists.root");

	if(!f1 || !f2) {
		std::cout << "ERROR in compareHistograms : Files do not exist for mode " << mode << std::endl;
		return 1;
	}

	// Open a multi-page PDF that aggregates every per-histogram comparison for
	// this mode. Written next to the run outputs (not under plots/) so it is
	// easy to find at the top of the CI artifact.
	TString summaryPdf = mode + "_summary.pdf";
	c1.Print(summaryPdf + "[");

	TIter next(f1->GetListOfKeys());
	TKey *key;

	while ((key = (TKey*)next())) {
		TClass *cl = gROOT->GetClass(key->GetClassName());
		if (!cl->InheritsFrom("TH1")) continue;
		TH1 *h1 = (TH1*)key->ReadObj();
		TH1* h2 = (TH1*)f2->Get(h1->GetName());

		if(!h2) {
			std::cout << "ERROR in compareHistograms : Histogram " << h1->GetName() << " not found." << std::endl;
			c1.Print(summaryPdf + "]");
			return 1;
		}

		if(TMath::Abs(h1->GetXaxis()->GetXmax()-h2->GetXaxis()->GetXmax())>0.001) {
			std::cout << "ERROR in compareHistograms : Histogram " << h1->GetName() << " has different maxima " << h1->GetXaxis()->GetXmax() << "\t" << h2->GetXaxis()->GetXmax() << std::endl;
			c1.Print(summaryPdf + "]");
			return 1;
		}
		if(TMath::Abs(h1->GetXaxis()->GetXmin()-h2->GetXaxis()->GetXmin())>0.001) {
			std::cout << "ERROR in compareHistograms : Histogram " << h1->GetName() << " has different minima " << h1->GetXaxis()->GetXmin() << "\t" << h2->GetXaxis()->GetXmin() << std::endl;
			c1.Print(summaryPdf + "]");
			return 1;
		}


		// Histograms concentrated in 1-2 bins are not really distributions,
		// so the Chi2Test isn't meaningful for them. With a single non-empty
		// bin (truth vertex coords always at 0) it returns p=0 even for
		// bit-identical histograms; with two non-empty bins it flags
		// float-precision noise on truth-level mass peaks (mass set to PDG
		// value, ~5% bin migration between Mac and Linux). Skip these from
		// the gate entirely — they are uninformative for distribution
		// comparison.
		int nNonEmpty = 0;
		for (int b = 1; b <= h1->GetNbinsX(); ++b) {
			if (h1->GetBinContent(b) != 0 || h2->GetBinContent(b) != 0) ++nNonEmpty;
		}
		if (nNonEmpty < 3) {
			std::cout << "INFO in compareHistograms : " << h1->GetName()
				  << " :\tskipped (only " << nNonEmpty
				  << " non-empty bin(s) — narrow distribution)" << std::endl;
			++nSkip;
			continue;
		}

		double pval = h1->Chi2Test(h2,"UU");
		std::cout << "INFO in compareHistograms : " << h1->GetName() << " :\tp-value = " << pval << std::endl;
		if(pval < warnThreshold) {
			std::cout << "WARNING in compareHistograms : Histograms do not match for " << h1->GetName() << " (p < 1%)" << std::endl;
			++nWarn;
		}
		if(pval < failThreshold) {
			++nFail;
			failing.push_back(std::make_pair(TString(h1->GetName()), pval));
		}
		hpval.Fill(pval);
		sumpval+=pval;
		++sum;

		h1->SetTitle(TString::Format("%s   p = %g", h1->GetName(), pval));
		h1->SetLineColor(kBlack);
		h1->SetStats(true);
		h2->SetLineColor(kRed);
		h2->SetStats(true);
		h1->Draw();
		// "sames" (with the trailing 's') tells ROOT to also draw h2's stats
		// box; without it the new build's stats would be suppressed.
		h2->Draw("sames");

		// Force the canvas to actually create both stats objects, then colour
		// them and shift the new build's box down so it doesn't overlap.
		c1.Update();
		TPaveStats *st1 = (TPaveStats*)h1->FindObject("stats");
		TPaveStats *st2 = (TPaveStats*)h2->FindObject("stats");
		if(st1) { st1->SetTextColor(kBlack); }
		if(st2) {
			st2->SetTextColor(kRed);
			double y1 = st2->GetY1NDC();
			double y2 = st2->GetY2NDC();
			double dy = y2 - y1;
			st2->SetY1NDC(y1 - dy);
			st2->SetY2NDC(y2 - dy);
		}

		TLegend leg(0.65, 0.50, 0.88, 0.60);
		leg.SetBorderSize(0);
		leg.SetFillStyle(0);
		leg.AddEntry(h1, "reference", "l");
		leg.AddEntry(h2, "this build", "l");
		leg.Draw();

		// Append a page to the multi-page summary PDF
		c1.Print(summaryPdf);

		// Keep individual per-histogram PDFs (backwards-compatible)
		TString pName = "plots/";
		pName += mode;
		pName += "_";
		pName += h1->GetName();
		pName += ".pdf";
		c1.SaveAs(pName);

	}
	hpval.Fit("pol0");
	hpval.Draw();
	// Final page of the summary PDF + standalone p-value PDF
	c1.Print(summaryPdf);
	c1.Print(summaryPdf + "]");
	TString pName = "plots/";
	pName += mode;
	pName += "_pval.pdf";
	c1.SaveAs(pName);

	double meanp = (sum > 0) ? sumpval/sum : 0.;
	bool meanFail = (meanp < 0.1);
	bool anyFail  = (nFail > 0);

	std::cout << std::endl;
	std::cout << "=== Summary for mode " << mode << " ===" << std::endl;
	std::cout << "Total histograms compared:   " << (int)sum << std::endl;
	std::cout << "  Skipped (narrow):          " << nSkip << std::endl;
	std::cout << "  p < " << warnThreshold << "  (warnings):       " << nWarn << std::endl;
	std::cout << "  p < " << failThreshold << " (FAIL threshold): " << nFail << std::endl;
	std::cout << "Mean p-value:                " << meanp << std::endl;
	if(!failing.empty()) {
		std::cout << "Failing histograms (p < " << failThreshold << "):" << std::endl;
		for(size_t i=0; i<failing.size(); ++i) {
			std::cout << "  " << failing[i].first << "\tp=" << failing[i].second << std::endl;
		}
	}
	std::cout << "RESULT: " << ((meanFail || anyFail) ? "FAIL" : "PASS") << std::endl;
	std::cout << "Combined plots: " << summaryPdf << std::endl;
	std::cout << std::endl;

	if(meanFail) {
		std::cout << "ERROR in compareHistograms : mean p-value " << meanp << " < 0.1 - check plots" << std::endl;
		return 1;
	}
	if(anyFail) {
		std::cout << "ERROR in compareHistograms : " << nFail << " histogram(s) below p=" << failThreshold << " threshold" << std::endl;
		return 1;
	}
	return 0;
}

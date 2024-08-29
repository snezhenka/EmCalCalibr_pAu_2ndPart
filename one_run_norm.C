#include <iostream>
#include <cmath>
#include <TRandom3.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH3.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include <TAxis.h>
#include <dirent.h>

void one_run_norm(const TString &directory)
{
    DIR *dir;
    struct dirent *ent;
    
    const int num_sectors = 8;
        std::vector<std::vector<double>> means(num_sectors);  // means[i] for sector i
        std::vector<std::vector<double>> errors(num_sectors); // errors[i] for sector i
        std::vector<std::vector<double>> segment_numbers(num_sectors); // segment numbers for plotting
    
    if ((dir = opendir(directory.Data())) != nullptr) {
        // Iterate over all files in the directory
        while ((ent = readdir(dir)) != nullptr) {
            TString filename = ent->d_name;
       
            std::string a = ent->d_name;
            a.erase(0,3);
            std::string outfilename("/home/pc/Trains/EmCalCalib_pAu/1DHistsFgBg/out_");
            outfilename.append(a);
            TFile *in_runs = new TFile(outfilename.c_str(), "RECREATE");

            
            // Check if the file matches the pattern *_se.root
            if (filename.BeginsWith("se-") &&
                filename.Length() == 14 &&
                filename.EndsWith(".root")) {    //если файл кончается на что-то то пишем EndsWith
                TString fullPath = directory + "/" + filename;
                TFile *file = TFile::Open(fullPath, "read");
    
    
  
    
    std::string hz("se-432639.root");
    hz.erase(0,3);
    hz.erase(6,5);
    hz.append("_sec_");

    std::string segment_number_str(filename(3, 6).Data()); // Extract segment number from filename
    int segment_number = std::stoi(segment_number_str);

    	
    	for (int sector = 0; sector<8; sector++)
    	{
            
            TH3D *bg = (TH3D*)file->Get("background")->Clone();
            TH3D *fg = (TH3D*)file->Get("foreground")->Clone();
    	int sec_bin = bg->GetYaxis()->FindBin((2*sector+1)/2);
    	
    	double sr;//srznach
   	double sigma;//sigma
            double error; //неопределенность
            
            in_runs->cd();
        
            std::string sektor = std::to_string(sector);
    	std::string histname_bg("bg_run_sec_"); //гистограмма одномерная для bg
    	histname_bg.append(sektor);
    	
    	
    	std::string histname_fg("fg_run_sec_"); //гистограмма одномерная для fg
    	histname_fg.append(sektor);
    	
    	TH1D* m_bg_hi_pt = (TH1D*)bg->ProjectionX(histname_bg.c_str(),sec_bin,sec_bin,2,10); //создание одномерной сразу из трехмерной,от 2 до 10 это импульс
    	TH1D* m_fg_hi_pt = (TH1D*)fg->ProjectionX(histname_fg.c_str(),sec_bin,sec_bin,2,10);
    


	TH1F *back = (TH1F*)m_bg_hi_pt->Clone();
	TH1F *razn = (TH1F*)m_fg_hi_pt->Clone();
	int left = m_fg_hi_pt->FindBin(0.1);
	int right = m_fg_hi_pt->FindBin(0.2);
	
	double int_fg_l = m_fg_hi_pt->Integral(1,left);
	double int_fg_r = m_fg_hi_pt->Integral(right,m_fg_hi_pt->GetNbinsX());
	double int_fg = int_fg_l + int_fg_r;
	
	double int_bg_l = m_bg_hi_pt->Integral(1,left);
	double int_bg_r = m_bg_hi_pt->Integral(right,m_bg_hi_pt->GetNbinsX());
	double int_bg = int_bg_l + int_bg_r;
	
	back->Scale(int_fg/int_bg);
	razn->Add(back,-1);
    	
    	hz.erase(11,1);
    	hz.append(sektor);
    	
    	 std::string fist_name("fist_run_sec_");
   	 fist_name.append(sektor);
   	 
    	 
   	 
   	 TF1 *fist_M = new TF1(hz.c_str(),"gaus(0)+pol2(3)", 0.05, 0.25);
   	 
   	 fist_M->SetParameter(0,1e6);
   	 fist_M->SetParameter(1,0.135);
   	 fist_M->SetParameter(2,0.015);
   	 
   	
   	 
    	 razn->Fit(fist_M, "0","QR");
    	 
            
    	 fist_M->SetParLimits(0,1,300e6);
   	 fist_M->SetParLimits(1,0.12,0.15);
   	 fist_M->SetParLimits(2,0.005,0.03);
   	 
   	 razn->Fit(fist_M, "QR");
    	 
    	 sr = fist_M -> GetParameter(1);//srznach
    	 error = fist_M -> GetParError(1); //неопределенность
   	 sigma = fist_M -> GetParameter(2);//sigma
      //  std::cout << "Mean: " << sr << "Error: " << error << std::endl;
            // Store the results for plotting
    means[sector].push_back(sr);
    errors[sector].push_back(error);
    segment_numbers[sector].push_back(segment_number);

            
            razn->Write();
        }
            file->Close();
            
         
        }
            in_runs->Close();
        }
    
                closedir(dir);
        
        // Plotting
               TFile *outFile = new TFile("/home/pc/Trains/EmCalCalib_pAu/1DHistsFgBg/graphs.root", "RECREATE");
               for (int sector = 0; sector < num_sectors; sector++) {
                   std::vector<Double_t> filtered_segment_numbers;
                       std::vector<Double_t> filtered_means;
                       std::vector<Double_t> filtered_errors;

                       for (size_t i = 0; i < segment_numbers[sector].size(); i++) {
                         //  if (errors[sector][i] <= 0.002 & errors[sector][i] > 0.00002 & means[sector][i] < 0.14 & means[sector][i] > 0.13) { // выставляем ограничение по ошибке и среднему
                               filtered_segment_numbers.push_back(segment_numbers[sector][i]);
                               filtered_means.push_back(means[sector][i]);
                               filtered_errors.push_back(errors[sector][i]);
                         //  }
                       }

                       if (filtered_segment_numbers.size() > 0) { //есть ли элементы
                           // Создание вектора пар data, первая часть сегменты данных,потом средние и их ошибки
                                  std::vector<std::pair<Double_t, std::pair<Double_t, Double_t>>> data;
                                  for (size_t i = 0; i < filtered_segment_numbers.size(); i++) {
                                      //заполнение вектора data
                                      data.push_back({filtered_segment_numbers[i], {filtered_means[i], filtered_errors[i]}});
                                  }

                                  // вектор data сортируется по первому элементу пар (по значениям сегмента данных),для сортировки исп лямбда функция которая сравнивает 2 элемента. нужно чтоб на графике точки соединялись последовательно
                                  std::sort(data.begin(), data.end(),
                                            [](const std::pair<Double_t, std::pair<Double_t, Double_t>> &a,
                                               const std::pair<Double_t, std::pair<Double_t, Double_t>> &b) {
                                                return a.first < b.first;
                                            });

                                  // после сортировки векторы очищаются и заполняются отсортированными значениями data
                                  filtered_segment_numbers.clear();
                                  filtered_means.clear();
                                  filtered_errors.clear();
                                  for (const auto &entry : data) {
                                      filtered_segment_numbers.push_back(entry.first);
                                      filtered_means.push_back(entry.second.first);
                                      filtered_errors.push_back(entry.second.second);
                                  }

                          
                           TGraphErrors *graph = new TGraphErrors(filtered_segment_numbers.size(),
                                                                        filtered_segment_numbers.data(),
                                                                        filtered_means.data(),
                                                                        0,
                                                                        filtered_errors.data());
                           graph->SetLineStyle(0);
                 /*  TGraphErrors *graph = new TGraphErrors(segment_numbers[sector].size(),
                                                           segment_numbers[sector].data(),
                                                           means[sector].data(),
                                                           0,
                                                           errors[sector].data()); */
                   graph->SetTitle(Form("Sector %d", sector));
                   graph->GetXaxis()->SetTitle("Segment Number");
                   graph->GetYaxis()->SetTitle("Mean Value");
                           graph->GetYaxis()->SetRangeUser(0.12,0.15);
                   // Fit a constant to the graph
                              TF1 *fit_const = new TF1(Form("const_fit_sector_%d", sector), "pol0", graph->GetXaxis()->GetXmin(), graph->GetXaxis()->GetXmax());
                           fit_const->SetParameter(0,0.135);
                           fit_const->SetParLimits(0,0.13,0.14);
                              graph->Fit(fit_const, "Q");
                           double constant = fit_const->GetParameter(0);
                           std::cout << "Constant of the fit: " << constant << std::endl;
                              fit_const->SetLineColor(kRed);
                              fit_const->SetLineWidth(2);
                           graph->SetLineWidth(1.5);
                          graph->SetMarkerSize(0.7);
                           graph->SetMarkerStyle(20); //круг
                           TCanvas *c1 = new TCanvas(Form("const_fit_sector_%d", sector), "Graph with Constant Fit", 800, 600);
                              graph->Draw("APE");
                              fit_const->Draw("same");
    
                   graph->Write(Form("graph_sector_%d", sector));
                           
               }
               }
               outFile->Close();

        }

}

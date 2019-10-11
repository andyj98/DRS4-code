/*
Short code to integrate waveforms and obtain histogram
This version of the code plots the amplitude of each peak, rather than integrating each peak
*/

#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH1.h"

using namespace std;

typedef struct {
   char           tag[3];
   char           version;
} FHEADER;

typedef struct {
   char           time_header[4];
} THEADER;

typedef struct {
   char           bn[2];
   unsigned short board_serial_number;
} BHEADER;

typedef struct {
   char           event_header[4];
   unsigned int   event_serial_number;
   unsigned short year;
   unsigned short month;
   unsigned short day;
   unsigned short hour;
   unsigned short minute;
   unsigned short second;
   unsigned short millisecond;
   unsigned short range;
} EHEADER;

typedef struct {
   char           tc[2];
   unsigned short trigger_cell;
} TCHEADER;

typedef struct {
   char           c[1];
   char           cn[3];
} CHEADER;

/*-----------------------------------------------------------------------------*/
//defining graph function
TGraph* graph(double a[0][4][1024]){

	int i;
//	TCanvas *c1 = new TCanvas("c1", "Channel 1 plotted waveforms", 200, 10, 700, 500);
	 
	/*c1->SetFillColor(42);
	c1->SetGrid(); */

	double x[1024], y[1024];

         /* channel 1 */
         for(i=0;i<1024;i++){
	   x[i] = i;
	   y[i] = a[0][0][i] * 1000;
         }

	 TGraph* gr = new TGraph(1024,x,y);
	 gr->SetLineColor(2);
	 gr->SetLineWidth(1);
	 gr->SetMarkerColor(4);
	 gr->SetMarkerStyle(21);
	 gr->SetTitle("Plotted waveforms");
	 gr->GetXaxis()->SetTitle("Voltage Bins");
	 gr->GetYaxis()->SetTitle("Voltage");
	 gr->Draw("ACP");
/*
	 c1->Update();
	 c1->Modified();
	 c1->Draw();

	 TFile* chn_graphs = TFile::Open("RECREATE");
	 chn_graphs->cd();
	 c1->Write("mycan");
	 c1->Draw();
	 chn_graphs->Close();
	
*/
	return gr;
	}

//defining functions to find max/min values of an array.
double findHighest( double a[][4][1024])
{
	double high = a[0][0][0];

	for(int i = 1; i < 1024; i++){
		if (a[0][0][i] > high){
			high = a[0][0][i]; }
			}
	return high;
}

double findLowest( double c[][4][1024])
{
	double low = c[0][0][0];

	for(int i = 1; i < 1024; i++){
		if (c[0][0][i] < low){
			low = c[0][0][i]; }
			}
	return low;
}
int main(int argc, const char * argv[])
{
   FHEADER  fh;
   THEADER  th;
   BHEADER  bh;
   EHEADER  eh;
   TCHEADER tch;
   CHEADER  ch;

   unsigned int scaler;
   unsigned short voltage[1024];
   double waveform[16][4][1024], time[16][4][1024];
   float bin_width[16][4][1024];
   int i, j, b, chn, n, chn_index, n_boards;
   double t1, t2, dt, temp1=0.0, temp2=0.0, temp3=0.0, temp4=0.0;
 //  double max=0.0, min=0.0;
   char filename[256];
   std::vector<double> waveform_area1, waveform_area2, waveform_area3, waveform_area4, height_chn1;

   int ndt;
   double threshold, sumdt, sumdt2;
   

   if (argc > 1)
      strcpy(filename, argv[1]);
   else {
      printf("Usage: read_binary <filename>\n");
      return 0;
   }

   // open the binary waveform file
   FILE *f = fopen(filename, "rb");
   if (f == NULL) {
      printf("Cannot find file \'%s\'\n", filename);
      return 0;
   }

   // read file header
   fread(&fh, sizeof(fh), 1, f);
   if (fh.tag[0] != 'D' || fh.tag[1] != 'R' || fh.tag[2] != 'S') {
      printf("Found invalid file header in file \'%s\', aborting.\n", filename);
      return 0;
   }

   if (fh.version != '2') {
      printf("Found invalid file version \'%c\' in file \'%s\', should be \'2\', aborting.\n", fh.version, filename);
      return 0;
   }

   // read time header
   fread(&th, sizeof(th), 1, f);
   if (memcmp(th.time_header, "TIME", 4) != 0) {
      printf("Invalid time header in file \'%s\', aborting.\n", filename);
      return 0;
   }

   for (b = 0 ; ; b++) {
      // read board header
      fread(&bh, sizeof(bh), 1, f);
      if (memcmp(bh.bn, "B#", 2) != 0) {
         // probably event header found
         fseek(f, -4, SEEK_CUR);
         break;
      }

      printf("Found data for board #%d\n", bh.board_serial_number);

      // read time bin widths
      memset(bin_width[b], sizeof(bin_width[0]), 0);
      for (chn=0 ; chn<5 ; chn++) {
         fread(&ch, sizeof(ch), 1, f);
         if (ch.c[0] != 'C') {
            // event header found
            fseek(f, -4, SEEK_CUR);
            break;
         }
         i = ch.cn[2] - '0' - 1;
         printf("Found timing calibration for channel #%d\n", i+1);
         fread(&bin_width[b][i][0], sizeof(float), 1024, f);
         // fix for 2048 bin mode: double channel
         if (bin_width[b][i][1023] > 10 || bin_width[b][i][1023] < 0.01) {
            for (j=0 ; j<512 ; j++)
               bin_width[b][i][j+512] = bin_width[b][i][j];
         }
      }
   }
   n_boards = b;



	string rootfile1;
        size_t pos2;
        rootfile1=argv[1];
        pos2=rootfile1.find('.');
        if(pos2!=rootfile1.npos){
          rootfile1=rootfile1.substr(0,pos2);
          rootfile1.append("_graph.root");
          cout<<"name of output graph file is "<<rootfile1<<endl;
        }

	const char* r1 = rootfile1.c_str();
//creating graphs and file to store saved waveforms
TGraph* gr = new TGraph(1024);
TFile* chn_graphs = TFile::Open(r1,"RECREATE");
	 chn_graphs->cd();


   // loop over all events in the data file
   for (n=0 ; ; n++) {
      // read event header
      i = (int)fread(&eh, sizeof(eh), 1, f);
      if (i < 1)
         break;

   //   printf("Found event #%d %d %d\n", eh.event_serial_number, eh.second, eh.millisecond);
	

      // loop over all boards in data file
      for (b=0 ; b<n_boards ; b++) {

         // read board header
         fread(&bh, sizeof(bh), 1, f);
         if (memcmp(bh.bn, "B#", 2) != 0) {
            printf("Invalid board header in file \'%s\', aborting.\n", filename);
            return 0;
         }

         // read trigger cell
         fread(&tch, sizeof(tch), 1, f);
         if (memcmp(tch.tc, "T#", 2) != 0) {
            printf("Invalid trigger cell header in file \'%s\', aborting.\n", filename);
            return 0;
         }

         if (n_boards > 1)
            printf("Found data for board #%d\n", bh.board_serial_number);

         // reach channel data
         for (chn=0 ; chn<4 ; chn++) {

            // read channel header
            fread(&ch, sizeof(ch), 1, f);
            if (ch.c[0] != 'C') {
               // event header found
               fseek(f, -4, SEEK_CUR);
               break;
            }
            chn_index = ch.cn[2] - '0' - 1;
            fread(&scaler, sizeof(int), 1, f);
            fread(voltage, sizeof(short), 1024, f);

            for (i=0 ; i<1024 ; i++) {
               // convert data to volts
               waveform[b][chn_index][i] = (voltage[i] / 65536. + eh.range/1000.0 - 0.5);

               // calculate time for this cell
               for (j=0,time[b][chn_index][i]=0 ; j<i ; j++)
                  time[b][chn_index][i] += bin_width[b][chn_index][(j+tch.trigger_cell) % 1024];
            }
         }


/* need to write each graph with each own specific name to prevent writing over existing objects */
	 //graphs all event waveforms 
	gr = graph(waveform);
	//gr->Draw();
	char name[15];
	sprintf(name, "event_%d", eh.event_serial_number);
	gr->Write(name);

/*Finding the height of each peak, and storing it into an array to be plotted instead of area under curve
	TGraph* max = TGraph::GetMaximum(name);
	TGraph* min = TGraph::GetMinimum(name);
	double diff = max - min;
	cout<< diff << endl; */


//Plotting the amplitude of each peak on chn. 1, rather than the integration of each one
	
         /* channel 1 */
         for(i=0;i<1024;i++){
           temp1+= waveform[b][0][i];
	//cout<<waveform[b][0][i]<<endl;
         }
	double max1 = findHighest(waveform);
	double min1 = findLowest(waveform);
	double diff1 = max1 - min1;
	height_chn1.push_back(diff1);
         waveform_area1.push_back(temp1);
         temp1=0;


         /* channel 2 */
         for(i=0;i<1024;i++){
           temp2+= waveform[b][1][i];
         }
          waveform_area2.push_back(temp2);
          temp2=0;

         /* channel 3 */
         for(i=0;i<1024;i++){
           temp3+= waveform[b][2][i];
         }
          waveform_area3.push_back(temp3);
          temp3=0;

         /* channel 4 */
         for(i=0;i<1024;i++){
           temp4+= waveform[b][3][i];
         }
         waveform_area4.push_back(temp4);
         temp4=0;

         // align cell #0 of all channels
         t1 = time[b][0][(1024-tch.trigger_cell) % 1024];
         for (chn=1 ; chn<4 ; chn++) {
            t2 = time[b][chn][(1024-tch.trigger_cell) % 1024];
            dt = t1 - t2;
            for (i=0 ; i<1024 ; i++)
               time[b][chn][i] += dt;
         }


	

         t1 = t2 = 0;
         threshold = 0.3;

         // find peak in channel 1 above threshold
         for (i=0 ; i<1022 ; i++)
            if (waveform[b][0][i] < threshold && waveform[b][0][i+1] >= threshold) {
               t1 = (threshold-waveform[b][0][i])/(waveform[b][0][i+1]-waveform[b][0][i])*(time[b][0][i+1]-time[b][0][i])+time[b][0][i];
               break;
            }

         // find peak in channel 2 above threshold
         for (i=0 ; i<1022 ; i++)
            if (waveform[b][1][i] < threshold && waveform[b][1][i+1] >= threshold) {
               t2 = (threshold-waveform[b][1][i])/(waveform[b][1][i+1]-waveform[b][1][i])*(time[b][1][i+1]-time[b][1][i])+time[b][1][i];
               break;
            }
          }
        }
//close graph file
chn_graphs->Close();

        TCanvas* can = new TCanvas("can","Raw Waveform Area",1200,1200);
      //  can->Divide(2,2,0.01,0.01);

        can->cd(1);
        //cout<<"***************** Channel 1 *****************"<<endl;
        TH1D* h1 = new TH1D("h1","Channel 1",1000,0,0.1);
        //unsigned int vector_size1 = waveform_area1.size();
	unsigned int vector_size1 = height_chn1.size();
        for(i=0;i<vector_size1;i++){
          //cout<<waveform_area1[i]<<endl;
          //h1->Fill(waveform_area1[i]);
	  h1->Fill(height_chn1[i]);
        }
	h1->SetTitle("Cs-137 Energy Spectrum");
	h1->GetXaxis()->SetTitle("Pulse height [V]");
	h1->GetYaxis()->SetTitle("Counts");
        h1->Draw();
        can->Modified();
        can->Update();

//Commenting out code to divide canvas into 4 histograms, one for each channel on the digitizer.//
 /*       can->cd(2);
        //cout<<"***************** Channel 2 *****************"<<endl;
        TH1D* h2 = new TH1D("h2","Channel 2",100,-50,50);
        unsigned int vector_size2 = waveform_area2.size();
        for(i=0;i<vector_size2;i++){
          //cout<<waveform_area2[i]<<endl;
          h2->Fill(waveform_area2[i]);
        }
	h2->SetTitle("Energy Spectrum: Channel 2");
	h2->GetXaxis()->SetTitle("Energy - arbitrary units");
	h2->GetYaxis()->SetTitle("Counts");
        h2->Draw("SAME");
        can->Modified();
        can->Update();

        can->cd(3);
        //cout<<"***************** Channel 3 *****************"<<endl;
        TH1D* h3 = new TH1D("h3","Channel 3",100,-50,50);
        unsigned int vector_size3 = waveform_area3.size();
        for(i=0;i<vector_size3;i++){
          //cout<<waveform_area3[i]<<endl;
          h3->Fill(waveform_area3[i]);
        }
	h3->SetTitle("Energy Spectrum: Channel 3");
	h3->GetXaxis()->SetTitle("Energy - arbitrary units");
	h3->GetYaxis()->SetTitle("Counts");
        h3->Draw("SAME");
        can->Modified();
        can->Update();

        can->cd(4);
        //cout<<"***************** Channel 4 *****************"<<endl;
        TH1D* h4 = new TH1D("h4","Channel 4",100,-50,50);
        unsigned int vector_size4 = waveform_area4.size();
        for(i=0;i<vector_size4;i++){
          //cout<<waveform_area4[i]<<endl;
          h4->Fill(waveform_area4[i]);
        }
	h4->SetTitle("Energy Spectrum: Channel 4");
	h4->GetXaxis()->SetTitle("Energy - arbitrary units");
	h4->GetYaxis()->SetTitle("Counts");
        h4->Draw("SAME");
        can->Modified();
        can->Update();
*/
        string rootfile;
        size_t pos1;
        rootfile=argv[1];
        pos1=rootfile.find('.');
        if(pos1!=rootfile.npos){
          rootfile=rootfile.substr(0,pos1);
          rootfile.append(".root");
          cout<<"name of output file is "<<rootfile<<endl;
        }

	

        const char* r = rootfile.c_str();
	

        TFile* fout = TFile::Open(r,"RECREATE");
        fout->cd();
        can->Write("can");
        can->Draw();
        fout->Close();
 
        return 0;
      }


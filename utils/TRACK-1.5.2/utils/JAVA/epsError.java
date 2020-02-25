package tigge;
import java.util.*;
import java.math.*;
import java.io.*;
import EPSAnalysis.*;

public class epsError{

        private static final int NUM_TRACKS_MAX = 1000000;

	public static void main(String[] args)throws IOException{
		
		
		String date = null;
		String filename = null;
		String fptmp = null;
		
		String[] fpAnalysis = new String[NUM_TRACKS_MAX];
		
		double[] tmpSepDiff = new double[Constants.LT_MAX];
		double[] tmpIntDiff = new double[Constants.LT_MAX];
		double[] tmpIntBias = new double[Constants.LT_MAX];
		double[] sumSepDiff = new double[Constants.LT_MAX];
		double[] sumIntDiff = new double[Constants.LT_MAX];
		double[] sumIntBias = new double[Constants.LT_MAX];
		double[] sumMeanSepDiff = new double[Constants.LT_MAX];
		double[] sumMeanIntDiff = new double[Constants.LT_MAX];
		double[] sumMeanIntBias = new double[Constants.LT_MAX];
		double[] sumCntlSepDiff = new double[Constants.LT_MAX];
		double[] sumCntlIntDiff = new double[Constants.LT_MAX];
		double[] sumCntlIntBias = new double[Constants.LT_MAX];
		
		double[] sumSepDiff2 = new double[Constants.LT_MAX];
		double[] sumIntDiff2 = new double[Constants.LT_MAX];
		double[] sumIntBias2 = new double[Constants.LT_MAX];
		double[] sumMeanSepDiff2 = new double[Constants.LT_MAX];
		double[] sumMeanIntDiff2 = new double[Constants.LT_MAX];
		double[] sumMeanIntBias2 = new double[Constants.LT_MAX];
		double[] sumCntlSepDiff2 = new double[Constants.LT_MAX];
		double[] sumCntlIntDiff2 = new double[Constants.LT_MAX];
		double[] sumCntlIntBias2 = new double[Constants.LT_MAX];
				
		int[] countMeanDiff = new int[Constants.LT_MAX];
		int[] countCntlDiff = new int[Constants.LT_MAX];
		int[] countDiff = new int[Constants.LT_MAX];
		
		TrackDate forSt;
		File ferr=null;
		int noOfTracks;
		int trackId;
		
		int ntr=0, icur=0;
		
		boolean test = false;
		
		BufferedWriter writer;
		
		PrintWriter[] cntrl_err = new PrintWriter[NUM_TRACKS_MAX];
		PrintWriter[] mean_err = new PrintWriter[NUM_TRACKS_MAX];
		PrintWriter[] eps_err = new PrintWriter[NUM_TRACKS_MAX];
		
// Read region file if present

                Region reg = RegionOp.ReadReg();
		
		Arrays.fill(cntrl_err, null);
		Arrays.fill(mean_err, null);
		Arrays.fill(eps_err, null);
		
                if(args.length != 7) {
                  System.err.println("Usage: epsError [path to dirs] [longitude id] [latitude id] [intensity id] [Match dir] [Mean Out Typ. 0 or 1] [Write Err., 0 or 1]");
                  System.exit(1);
                }

               
		int ilong = Integer.parseInt(args[1]);
		int ilat = Integer.parseInt(args[2]);
		int ifd = Integer.parseInt(args[3]);
		int ipt = Integer.parseInt(args[5]);
		int iwe = Integer.parseInt(args[6]);
		
		if(ipt < 0 || ipt > 1) {
		   System.err.println("Mean output type specifier incorrect, must be 0 or 1");
		   System.exit(1);
		}
		if(iwe < 0 || iwe > 1) {
		   System.err.println("Mean output type specifier incorrect, must be 0 or 1");
		   System.exit(1);
		}		
		

                String matchd = args[4];
		
                File[] forDates = (new File(args[0])).listFiles();
		Arrays.sort(forDates);

                System.out.println(forDates.length);
		
		FilenameFilter fmatchtyp = new FilenameFilter() {  
                   public boolean accept(File d, String name) { 
		      boolean ret;
		      if(name != null && (name.startsWith("trmatch") && ! name.endsWith("_mean"))) {
		        ret = true;
		      }
		      else ret = false;
		      return ret;
		   }
		};
		
		
		for(int i=0; i<Constants.LT_MAX; i++){
				
			tmpSepDiff[i]=-100;
			tmpIntDiff[i]=-100;
			tmpIntBias[i]=-100;
				
			sumMeanSepDiff[i]=0;
			sumMeanIntBias[i]=0;
			sumMeanIntDiff[i]=0;

			sumMeanSepDiff2[i]=0;
			sumMeanIntBias2[i]=0;
			sumMeanIntDiff2[i]=0;			
							
			sumCntlSepDiff[i]=0;
			sumCntlIntBias[i]=0;
			sumCntlIntDiff[i]=0;
			
			sumCntlSepDiff2[i]=0;
			sumCntlIntBias2[i]=0;
			sumCntlIntDiff2[i]=0;			
				
			sumSepDiff[i]=0;
			sumIntBias[i]=0;								
			sumIntDiff[i]=0;
			
			sumSepDiff2[i]=0;
			sumIntBias2[i]=0;								
			sumIntDiff2[i]=0;			

                        countMeanDiff[i]=0;
			countCntlDiff[i]=0;
			countDiff[i]=0;
				
		}
		
		for(int i=0; i<forDates.length; i++){
		
			System.out.println(i + " of " + forDates.length);
			System.out.println(forDates[i].toString());
							
			date = forDates[i].toString();
                        date = date.substring(date.length() - 10);

			forSt = new TrackDate(Integer.parseInt(date.substring(0,4)), Integer.parseInt(date.substring(4,6)), Integer.parseInt(date.substring(6,8)), Integer.parseInt(date.substring(8,10)));			

			File curdir = new File(forDates[i].toString() + "/" + matchd + "/");
			
			//System.out.println(curdir.toString());
			
			if (curdir.isDirectory()) {
			
		            File[] trackFile = curdir.listFiles(fmatchtyp);
			    Arrays.sort(trackFile);
			    
                            if(trackFile != null){
			       //System.out.println(trackFile.length);
			 		
			       for(int j=0; j<trackFile.length; j++){
			          for(int l=0; l<Constants.LT_MAX; l++){
				     tmpSepDiff[l]=-100;
				     tmpIntDiff[l]=-100;
				     tmpIntBias[l]=-100;
				  }
				  
				  //System.out.println(trackFile[j]);
				 
// read analysis file

				  Track analysis = TrackOperations.File2Track(trackFile[j], 0, ilong, ilat, ifd);//analysis track
				  						
				  noOfTracks = TrackOperations.getNoOfTracks(trackFile[j]);
				  
// apply region masking to verification track

                                  if(reg.getIreg() == 1) analysis = TrackOperations.TrackMaskRegion(analysis, reg);				  
				  
// allocate analysis identifiers

                                  if(iwe == 1){

                                     fptmp = analysis.getPoint(0).getTrackDate() + "\t" + analysis.getPoint(0).getLon() + "\t" + analysis.getPoint(0).getLat();
				     //System.out.println(fptmp);
				  	
                                     if(!Arrays.asList(fpAnalysis).contains(fptmp)){
				        fpAnalysis[ntr] = fptmp;
				        ++ntr;		      
				     }
				  
				     icur = Arrays.asList(fpAnalysis).indexOf(fptmp);				  
				     //System.out.println(fptmp + "\t" + fpAnalysis[icur]);
				  
// open file for writing track error data
				  
				  
				     if(eps_err[icur] == null) {
				        ferr = new File("TRERROR/eps_err" + icur + ".dat");
				        eps_err[icur] = new PrintWriter(new BufferedWriter(new FileWriter(ferr, true)), true);
				     }
				  
				  }
				  
// read mean track file

                                  File meanfile = new File(trackFile[j].toString() + "_mean");

                                  if(meanfile.exists()){ 				  

                                     Track mean = TrackOperations.File2Track(meanfile,0, ilong, ilat, ifd);//mean track
				     tmpSepDiff = TrackOperations.separationDistance(analysis,mean,forSt);				  
				     tmpIntDiff = TrackOperations.intensityDifference(analysis,mean,forSt);						
				     tmpIntBias = TrackOperations.intensityDifference2(analysis,mean,forSt);
				     
				     if(iwe == 1){				  
                                        if(mean_err[icur] == null) {
				           ferr = new File("TRERROR/mean_err" + icur + ".dat");
				           mean_err[icur] = new PrintWriter(new BufferedWriter(new FileWriter(ferr, true)), true);
				        }
					
					mean_err[icur].printf("TTT %s %d %n", forSt, 0);

				     }
				  
// compute mean seperation difference and intensity differences for the mean track
				  
				     test = false;
				     for(int l=0; l<Constants.LT_MAX; l++){
							
				        if(tmpSepDiff[l] != -100 && tmpIntDiff[l] != -100 && tmpIntBias[l] != -100){
					   sumMeanSepDiff[l]+=tmpSepDiff[l];
					   sumMeanIntDiff[l]+=tmpIntDiff[l];					
					   sumMeanIntBias[l]+=tmpIntBias[l];
					
					   sumMeanSepDiff2[l]+=tmpSepDiff[l] * tmpSepDiff[l];
					   sumMeanIntDiff2[l]+=tmpIntDiff[l] * tmpIntDiff[l];					
					   sumMeanIntBias2[l]+=tmpIntBias[l] * tmpIntBias[l];
					   
					   if(iwe == 1) mean_err[icur].printf("%4d %e %e %e %n", l, tmpSepDiff[l], tmpIntDiff[l], tmpIntBias[l]);				   
										
					   countMeanDiff[l]++;
					   test = true;
		
				        }

				     }//end of l for loop

                                  }
				  
// determine if a control forecast is present and compute seperation and intensity differences.

                                  trackId = TrackOperations.getTrackID(trackFile[j], 1);
				  //System.out.println("TRACK_ID=" + "\t" + trackId);
                                  if(trackId == 1){

                                     for(int l=0; l<Constants.LT_MAX; l++){
                                        tmpSepDiff[l]=-100;
                                        tmpIntDiff[l]=-100;
                                        tmpIntBias[l]=-100;
                                     }
				     
				     Track control = TrackOperations.File2Track(trackFile[j], 1, ilong, ilat, ifd);//control track			     
				     tmpSepDiff = TrackOperations.separationDistance(analysis,control,forSt);
				     tmpIntDiff = TrackOperations.intensityDifference(analysis,control,forSt);						
				     tmpIntBias = TrackOperations.intensityDifference2(analysis,control,forSt);
				     
				     if(iwe == 1){				  
                                        if(cntrl_err[icur] == null) {
				           ferr = new File("TRERROR/cntrl_err" + icur + ".dat");
				           cntrl_err[icur] = new PrintWriter(new BufferedWriter(new FileWriter(ferr, true)), true);
				        }
					
					cntrl_err[icur].printf("TTT %s %d %n", forSt, 0);

				     }
				     
				     test = false;
				     for(int l=0; l<Constants.LT_MAX; l++){
							
				        if(tmpSepDiff[l] != -100 && tmpIntDiff[l] != -100 && tmpIntBias[l] != -100){
				   	   sumCntlSepDiff[l]+=tmpSepDiff[l];
					   sumCntlIntDiff[l]+=tmpIntDiff[l];					
					   sumCntlIntBias[l]+=tmpIntBias[l];
					   
				   	   sumCntlSepDiff2[l]+=tmpSepDiff[l] * tmpSepDiff[l];
					   sumCntlIntDiff2[l]+=tmpIntDiff[l] * tmpIntDiff[l];					
					   sumCntlIntBias2[l]+=tmpIntBias[l] * tmpIntBias[l];
					   
					   if(iwe == 1) cntrl_err[icur].printf("%4d %e %e %e %n", l, tmpSepDiff[l], tmpIntDiff[l], tmpIntBias[l]);
					   					   
					   countCntlDiff[l]++;
					   test = true;

							
				        }

				     }
				  
				  }
				  
// compute seperation and intensity differences over all ensemble members


                                  for(int k=0; k < noOfTracks - 1; k++){
                                     for(int l=0; l<Constants.LT_MAX; l++){
                                         tmpSepDiff[l]=-100;
                                         tmpIntDiff[l]=-100;
                                         tmpIntBias[l]=-100;
                                      }

                                      if(k == 0){
                                        trackId = TrackOperations.getTrackID(trackFile[j], 1);
                                        if(trackId == 1) continue;
                                      }

				      Track enstr = TrackOperations.File2Track(trackFile[j], k + 1, ilong, ilat, ifd);//ensemble track excluding control
				      tmpSepDiff = TrackOperations.separationDistance(analysis,enstr,forSt);
				      tmpIntDiff = TrackOperations.intensityDifference(analysis,enstr,forSt);						
				      tmpIntBias = TrackOperations.intensityDifference2(analysis,enstr,forSt);
				      
				      if(iwe == 1){
				         trackId = TrackOperations.getTrackID(trackFile[j], k + 1);
				         eps_err[icur].printf("TTT %s %d %n", forSt, trackId - 2);
				      }
				      
				      test = false;
				      for(int l=0; l<Constants.LT_MAX; l++){
				      
				        if(tmpSepDiff[l] != -100 && tmpIntDiff[l] != -100 && tmpIntBias[l] != -100){
				   	   sumSepDiff[l]+=tmpSepDiff[l];
					   sumIntDiff[l]+=tmpIntDiff[l];					
					   sumIntBias[l]+=tmpIntBias[l];
					   
				   	   sumSepDiff2[l]+=tmpSepDiff[l] * tmpSepDiff[l];
					   sumIntDiff2[l]+=tmpIntDiff[l] * tmpIntDiff[l];					
					   sumIntBias2[l]+=tmpIntBias[l] * tmpIntBias[l];
					   
					   if(iwe == 1) eps_err[icur].printf("%4d %e %e %e %n", l, tmpSepDiff[l], tmpIntDiff[l], tmpIntBias[l]);
					   					   
					   countDiff[l]++;
					   test = true;
		
				        }

				     }			      
				  
				  }		  	       
				  
			       }
			
			    }
			
			}
			
			
		}//end of i for loop
		
		
// close error files

                if(iwe == 1){
                   for(int i=0; i<NUM_TRACKS_MAX; i++){
		      if(cntrl_err[i] != null) cntrl_err[i].close();
		      if(mean_err[i] != null) mean_err[i].close();
		      if(eps_err[i] != null) eps_err[i].close();
		   }
		}
		
// output mean stats

		if(ipt == 0) writer = new BufferedWriter(WriteFile.openfile(null));
		else writer = new BufferedWriter(WriteFile.openfile("mean.txt")); 
		
		writer.write("Mean Track: Forecast points, Seperation Distance, Intensity Difference, Intensity Bias, Seperation STD, Intensity STD, Intensity Bias STD");
		writer.newLine();
		writer.flush();
					
		for(int i=0; i<Constants.LT_MAX; i++){
		        sumMeanSepDiff[i]/=countMeanDiff[i];
			sumMeanIntDiff[i]/=countMeanDiff[i];
			sumMeanIntBias[i]/=countMeanDiff[i];
		        sumMeanSepDiff2[i]/=countMeanDiff[i];
			sumMeanIntDiff2[i]/=countMeanDiff[i];
			sumMeanIntBias2[i]/=countMeanDiff[i];			
			writer.write(Round.round((double)i/4,2) + "\t" + countMeanDiff[i] + "\t" 
			                 + Round.round(sumMeanSepDiff[i],4) + "\t" 
					 + Round.round(sumMeanIntDiff[i],4) + "\t" 
					 + Round.round(sumMeanIntBias[i],4) + "\t"
					 + Round.round(Math.sqrt(sumMeanSepDiff2[i] - sumMeanSepDiff[i] * sumMeanSepDiff[i]),4) + "\t"
					 + Round.round(Math.sqrt(sumMeanIntDiff2[i] - sumMeanIntDiff[i] * sumMeanIntDiff[i]),4) + "\t"
					 + Round.round(Math.sqrt(sumMeanIntBias2[i] - sumMeanIntBias[i] * sumMeanIntBias[i]),4)        );
			writer.newLine();
		        writer.flush();
		}
		
		if(ipt == 1) writer.close();
		
// output control stats

		if(ipt == 0) writer = new BufferedWriter(WriteFile.openfile(null));
		else writer = new BufferedWriter(WriteFile.openfile("control.txt")); 

		writer.write("Control Track: Forecast points, Seperation Distance, Intensity Difference, Intensity Bias, Seperation STD, Intensity STD, Intensity Bias STD");
	        writer.newLine();
		writer.flush();
					
		for(int i=0; i<Constants.LT_MAX; i++){
		        sumCntlSepDiff[i]/=countCntlDiff[i];
			sumCntlIntDiff[i]/=countCntlDiff[i];
			sumCntlIntBias[i]/=countCntlDiff[i];
			sumCntlSepDiff2[i]/=countCntlDiff[i];
			sumCntlIntDiff2[i]/=countCntlDiff[i];
			sumCntlIntBias2[i]/=countCntlDiff[i];
			writer.write(Round.round((double)i/4,2) + "\t" + countCntlDiff[i] + "\t" 
			                 + Round.round(sumCntlSepDiff[i],4) + "\t" 
					 + Round.round(sumCntlIntDiff[i],4) + "\t" 
					 + Round.round(sumCntlIntBias[i],4) + "\t"
					 + Round.round(Math.sqrt(sumCntlSepDiff2[i] - sumCntlSepDiff[i] * sumCntlSepDiff[i]),4) + "\t"
					 + Round.round(Math.sqrt(sumCntlIntDiff2[i] - sumCntlIntDiff[i] * sumCntlIntDiff[i]),4) + "\t"
					 + Round.round(Math.sqrt(sumCntlIntBias2[i] - sumCntlIntBias[i] * sumCntlIntBias[i]),4)		);
			writer.newLine();
			writer.flush();
		}
		
		if(ipt == 1) writer.close();
		
		// output ensemble stats
		
		if(ipt == 0) writer = new BufferedWriter(WriteFile.openfile(null));
		else writer = new BufferedWriter(WriteFile.openfile("eps.txt"));
		
		writer.write("Ensemble Track: Forecast points, Seperation Distance, Intensity Difference, Intensity Bias, Seperation STD, Intensity STD, Intensity Bias STD");
	        writer.newLine();
		writer.flush();
					
		for(int i=0; i<Constants.LT_MAX; i++){
		        sumSepDiff[i]/=countDiff[i];
			sumIntDiff[i]/=countDiff[i];
			sumIntBias[i]/=countDiff[i];
		        sumSepDiff2[i]/=countDiff[i];
			sumIntDiff2[i]/=countDiff[i];
			sumIntBias2[i]/=countDiff[i];			
			writer.write(Round.round((double)i/4,2) + "\t" + countDiff[i] + "\t" 
			                 + Round.round(sumSepDiff[i],4) + "\t" 
					 + Round.round(sumIntDiff[i],4) + "\t" 
					 + Round.round(sumIntBias[i],4) + "\t"
					 + Round.round(Math.sqrt(sumSepDiff2[i] - sumSepDiff[i] * sumSepDiff[i]),4) + "\t"
					 + Round.round(Math.sqrt(sumIntDiff2[i] - sumIntDiff[i] * sumIntDiff[i]),4) + "\t"
					 + Round.round(Math.sqrt(sumIntBias2[i] - sumIntBias[i] * sumIntBias[i]),4)		 );
			writer.newLine();
			writer.flush();
		}
		
		if(ipt == 1) writer.close();
			
	}//end of main method
}//end of class

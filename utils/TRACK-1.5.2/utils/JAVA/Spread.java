package tigge;
import java.util.*;
import java.math.*;
import java.io.*;
import EPSAnalysis.*;

public class Spread{

        private static final int NUM_TRACKS_MAX = 1000000;
	
	public static void main(String[] args)throws IOException{
		
		String date = null;
		String filename = null;
		String fptmp = null;
		
		String[] fpAnalysis = new String[NUM_TRACKS_MAX];
		
		double[] tmpSepDiff = new double[Constants.LT_MAX];
		double[] tmpIntDiff = new double[Constants.LT_MAX];
		
		
		double[] tmpAnSepDiff = new double[Constants.LT_MAX];
		double[] tmpAnIntDiff = new double[Constants.LT_MAX];
				
		double[] sumSepDiff2 = new double[Constants.LT_MAX];
		double[] sumIntDiff2 = new double[Constants.LT_MAX];
		double[] spreadSep = new double[Constants.LT_MAX];
		double[] spreadInt = new double[Constants.LT_MAX];
		
		double[] spreadSep2 = new double[Constants.LT_MAX];
		double[] spreadInt2 = new double[Constants.LT_MAX];	
		
		double[] SepDiffMin = new double[Constants.LT_MAX];
		double[] SepDiffMax = new double[Constants.LT_MAX];
		
		double[] IntDiffMin = new double[Constants.LT_MAX];
		double[] IntDiffMax = new double[Constants.LT_MAX];	
		
		double[] sumSepMRE = new double[Constants.LT_MAX];
		double[] sumIntMRE = new double[Constants.LT_MAX];
		
		double[][] sumSepPWS = new double[3][Constants.LT_MAX];
		double[][] sumIntPWS = new double[3][Constants.LT_MAX];

		
		int[] countDiff2 = new int[Constants.LT_MAX];
		int[] spreadCount= new int[Constants.LT_MAX];
		int[] MREcount= new int[Constants.LT_MAX];
		int[] PWScount= new int[Constants.LT_MAX];		
		
		TrackDate forSt;
		File analFile;
		int noOfTracks;
		File ferr=null;
		
		double trat1, trat2;
		double MREcor;
		
		int ntr=0, icur=0;

		BufferedWriter writer;
		
		PrintWriter[] sprd_err = new PrintWriter[NUM_TRACKS_MAX];
		
// Read region file if present

                Region reg = RegionOp.ReadReg();
		
		Arrays.fill(sprd_err, null);
		
                if(args.length != 9) {
                  System.err.println("Usage: Spread [path to dirs] [longitude id] [latitude id] [intensity id] [Match dir] [Mean Out Typ. 0 or 1] [Write Err., 0 or 1] [Minimum No. Memb.] [Ensemble Size + control]");
                  System.exit(1);
                }

               
		int ilong = Integer.parseInt(args[1]);
		int ilat = Integer.parseInt(args[2]);
		int ifd = Integer.parseInt(args[3]);
		int ipt = Integer.parseInt(args[5]);
		int iwe = Integer.parseInt(args[6]);
                int nmmin = Integer.parseInt(args[7]); 
		int NeMax = Integer.parseInt(args[8]);

                if(nmmin <= 0) {
                   System.err.println("Minimum number of ensemble members must be greater than zero.");
                   System.exit(1);
                }

		
		NeMax /= 2;

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
				tmpAnSepDiff[i]=-100;
				tmpAnIntDiff[i]=-100;
				SepDiffMin[i]=-100;
				SepDiffMax[i]=-100;
				IntDiffMin[i]=-100;
				IntDiffMax[i]=-100;				
				sumSepDiff2[i]=0;
				sumIntDiff2[i]=0;
				sumSepMRE[i]=0;
				sumIntMRE[i]=0;
				sumSepPWS[0][i]=0;
				sumSepPWS[1][i]=0;
				sumSepPWS[2][i]=0;
				sumIntPWS[0][i]=0;
				sumIntPWS[1][i]=0;
				sumIntPWS[2][i]=0;				
				countDiff2[i]=0;
				spreadSep[i]=0;
				spreadInt[i]=0;
				spreadSep2[i]=0;
				spreadInt2[i]=0;				
				spreadCount[i]=0;
				MREcount[i]=0;
				PWScount[i]=0;
				
		}
		
		for(int i=0; i<forDates.length; i++){
		
			System.out.println(i + " of " + forDates.length);
						
			//System.out.println(forDates[i].toString());
			
			
			date = forDates[i].toString();
                        date = date.substring(date.length() - 10);
				
			//System.out.println(date);
				
			forSt = new TrackDate(Integer.parseInt(date.substring(0,4)), Integer.parseInt(date.substring(4,6)), Integer.parseInt(date.substring(6,8)), Integer.parseInt(date.substring(8,10)));

			//System.out.println(forSt);			
			
			File curdir = new File(forDates[i].toString() + "/" + matchd + "/");

                        //System.out.println(curdir.toString());

                        if (curdir.isDirectory()) {

                            File[] trackFile = curdir.listFiles(fmatchtyp);
			    Arrays.sort(trackFile);	
			    
			    if(trackFile != null){
		
				for(int j=0; j<trackFile.length; j++){
				
				    for(int l=0; l<Constants.LT_MAX; l++){
                                       tmpSepDiff[l]=-100;
                                       tmpIntDiff[l]=-100;

                                    }
				
					
// read mean track file

                                    File meanfile = new File(trackFile[j].toString() + "_mean");
				    
				    if(!meanfile.exists()) continue;

                                    Track mean = TrackOperations.File2Track(meanfile,0, ilong, ilat, ifd);//mean track
				    
// apply region masking to verification track

                                    //if(reg.getIreg() == 1) mean = TrackOperations.TrackMaskRegion(mean, reg);   

                                    //for(int l=0; l<mean.getPointNo(); l++){
                                    //   if(mean.getPoint(l).getIntensity() > 100){
                                    //      System.out.println(meanfile.toString());
                                    //   }

                                    //}
				    
 // read analysis file

                                    Track analysis = TrackOperations.File2Track(trackFile[j], 0, ilong, ilat, ifd);//analysis track
				    
				    noOfTracks = TrackOperations.getNoOfTracks(trackFile[j]);
				    
// apply region masking to verification track and mean track

                                   if(reg.getIreg() == 1) {
				   
				      analysis = TrackOperations.TrackMaskRegion(analysis, reg);
				      mean = TrackOperations.TrackMaskTime(analysis, mean);
				   }
				   				    
				    
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
				  
				  
				       if(sprd_err[icur] == null) {
				          ferr = new File("TRERROR/sprd_err" + icur + ".dat");
				          sprd_err[icur] = new PrintWriter(new BufferedWriter(new FileWriter(ferr, true)), true);
				       }
				       
				       sprd_err[icur].printf("TTT %s %d %n", forSt, 0);
				  
				    }
				    
				    
				    for(int l=0; l<Constants.LT_MAX; l++){
								
					tmpSepDiff[l]=-100;
					tmpIntDiff[l]=-100;
					tmpAnSepDiff[l]=-100;
					tmpAnIntDiff[l]=-100;
				        SepDiffMin[l]=-100;
				        SepDiffMax[l]=-100;
				        IntDiffMin[l]=-100;
				        IntDiffMax[l]=-100;					
					sumSepDiff2[l]=0;
					sumIntDiff2[l]=0;
					countDiff2[l]=0;
						
				    }
				    
				    tmpAnSepDiff = TrackOperations.separationDistance(mean,analysis,forSt);						
							
				    tmpAnIntDiff = TrackOperations.intensityDifference(mean,analysis,forSt);
						
//Calculate Spread
						
				    for(int l=1; l<noOfTracks; l++){ 
						
						
					Track ensMember = TrackOperations.File2Track(trackFile[j], l, ilong, ilat, ifd);//ensemble member track
							
					tmpSepDiff = TrackOperations.separationDistance(mean,ensMember,forSt);						
							
					tmpIntDiff = TrackOperations.intensityDifference(mean,ensMember,forSt);
							
					for(int k=0; k<Constants.LT_MAX; k++){
					    if(tmpSepDiff[k] != -100 && tmpIntDiff[k] != -100){			
						sumSepDiff2[k]+=tmpSepDiff[k];
						sumIntDiff2[k]+=tmpIntDiff[k];
						countDiff2[k]++;							
				  	    }
					  
// Calculate min/max range for MRE					  
					    if(SepDiffMin[k] == -100) SepDiffMin[k] = tmpSepDiff[k];
					    else if(tmpSepDiff[k] < SepDiffMin[k]) SepDiffMin[k] = tmpSepDiff[k];
					    
					    if(SepDiffMax[k] == -100) SepDiffMax[k] = tmpSepDiff[k];
					    else if(tmpSepDiff[k] > SepDiffMax[k]) SepDiffMax[k] = tmpSepDiff[k];

					    if(IntDiffMin[k] == -100) IntDiffMin[k] = tmpIntDiff[k];
					    else if(tmpIntDiff[k] < IntDiffMin[k]) IntDiffMin[k] = tmpIntDiff[k];
					    
					    if(IntDiffMax[k] == -100) IntDiffMax[k] = tmpIntDiff[k];
					    else if(tmpIntDiff[k] > IntDiffMax[k]) IntDiffMax[k] = tmpIntDiff[k];
					    					    
					    

					}//end of k loop
												
						
				     }//end  of l loop
				     
						
				     for(int l=0; l<Constants.LT_MAX; l++){
						
							
					if(countDiff2[l] >= nmmin){
					   trat1 = sumSepDiff2[l]/countDiff2[l];
					   trat2 = sumIntDiff2[l]/countDiff2[l];
					   spreadSep[l]+=trat1;
					   spreadInt[l]+=trat2;

					   spreadSep2[l]+=trat1 * trat1;
					   spreadInt2[l]+=trat2 * trat2;
					   
					   if(iwe == 1) sprd_err[icur].printf("%4d %e %e %e %n", l, trat1, trat2, 0.0);
							
					   spreadCount[l]++;

// Calculate MRE
					   if(countDiff2[l] >= NeMax){
					     if(SepDiffMin[l] != -100 && IntDiffMin[l] != -100) {
					       if(tmpAnSepDiff[l] != -100 && tmpAnIntDiff[l] != 100){
					         if( tmpAnSepDiff[l] < SepDiffMin[l] || tmpAnSepDiff[l] > SepDiffMax[l]) sumSepMRE[l] += 1.0;
					         if( tmpAnIntDiff[l] < IntDiffMin[l] || tmpAnIntDiff[l] > IntDiffMax[l]) sumIntMRE[l] += 1.0;
					         sumSepMRE[l] -= 2.0 / (countDiff2[l] + 1);
					         sumIntMRE[l] -= 2.0 / (countDiff2[l] + 1);
					         MREcount[l]++;
					       }
					     }
					     
					     if(tmpAnSepDiff[l] != -100 && tmpAnIntDiff[l] != 100){
					        if(tmpAnSepDiff[l] <= trat1) sumSepPWS[0][l] += 1.0;
						if(tmpAnSepDiff[l] <= 2.0 * trat1) sumSepPWS[1][l] += 1.0;
						if(tmpAnSepDiff[l] <= 3.0 * trat1) sumSepPWS[2][l] += 1.0;
						
					        if(tmpAnIntDiff[l] <= trat2) sumIntPWS[0][l] += 1.0;
						if(tmpAnIntDiff[l] <= 2.0 * trat2) sumIntPWS[1][l] += 1.0;
						if(tmpAnIntDiff[l] <= 3.0 * trat2) sumIntPWS[2][l] += 1.0;						
					     
					        PWScount[l]++;
					     
					     }
					     
					   }		
																								
					}
				     }//end of l for loop
						
				  
			       }
			       
			    }		
		
			}//end of j for loop
		
		}//end of i for loop
		
		
		for(int i=0; i<Constants.LT_MAX; i++){
//		   System.out.println(spreadCount[i] + "\t" + MREcount[i]);
                   if(MREcount[i] != 0){
		     sumSepMRE[i] = sumSepMRE[i] / MREcount[i];
		     sumIntMRE[i] = sumIntMRE[i] / MREcount[i];
		   }
		   if(PWScount[i] != 0){
		      sumSepPWS[0][i] /= PWScount[i];
		      sumSepPWS[1][i] /= PWScount[i];
		      sumSepPWS[2][i] /= PWScount[i];		   
		      sumIntPWS[0][i] /= PWScount[i];
		      sumIntPWS[1][i] /= PWScount[i];
		      sumIntPWS[2][i] /= PWScount[i];
		   }
		}
		
		writer = new BufferedWriter(WriteFile.openfile("MRE.txt"));
		writer.write("MRE: MRE Points, MRE Separation, MRE Intensity");
		writer.newLine();
		writer.flush();
		for(int i=0; i<Constants.LT_MAX; i++){		
		    writer.write(Round.round((double)i/4,2) + "\t" + MREcount[i] + "\t"
		                 + Round.round(sumSepMRE[i],4) + "\t"
				 + Round.round(sumIntMRE[i],4)                            );
		    writer.newLine();
		    writer.flush();
		}
		writer.close();
		
		writer = new BufferedWriter(WriteFile.openfile("PWS.txt"));
		writer.write("PWS: PWS Points, Separation1, Separation2, Separation3, Intensity1, Intensity2, Intensity3");
		writer.newLine();
		writer.flush();
		for(int i=0; i<Constants.LT_MAX; i++){		
		    writer.write(Round.round((double)i/4,2) + "\t" + PWScount[i] + "\t"
		                 + Round.round(sumSepPWS[0][i],4) + "\t" + Round.round(sumSepPWS[1][i],4) + "\t" + Round.round(sumSepPWS[2][i],4) + "\t"
				 + Round.round(sumIntPWS[0][i],4) + "\t" + Round.round(sumIntPWS[1][i],4) + "\t" + Round.round(sumIntPWS[2][i],4));
		    writer.newLine();
		    writer.flush();
		}
		writer.close();
			
		// close error files

                if(iwe == 1){
                   for(int i=0; i<NUM_TRACKS_MAX; i++){
		      if(sprd_err[i] != null) sprd_err[i].close();		      
		   }
		}
		
	        if(ipt == 0) writer = new BufferedWriter(WriteFile.openfile(null));
		else writer = new BufferedWriter(WriteFile.openfile("spread.txt")); 

                writer.write("Spread: Spread Points, Spread Separation, Spread Intensity, Seperation STD, Intensity STD");
		writer.newLine();
		writer.flush();
		
		for(int i=0; i<Constants.LT_MAX; i++){
		        spreadSep[i]/=spreadCount[i];
			spreadInt[i]/=spreadCount[i];
		        spreadSep2[i]/=spreadCount[i];
			spreadInt2[i]/=spreadCount[i];			
			writer.write(Round.round((double)i/4,2) + "\t" + spreadCount[i] + "\t" 
			                 + Round.round(spreadSep[i],4) + "\t" 
					 + Round.round(spreadInt[i],4) + "\t"
					 + Round.round(Math.sqrt(spreadSep2[i] - spreadSep[i] * spreadSep[i]),4) + "\t"
					 + Round.round(Math.sqrt(spreadInt2[i] - spreadInt[i] * spreadInt[i]),4)	 );
			writer.newLine();
		        writer.flush();
		}
		
		if(ipt == 1) writer.close();
		
		
	}//end of main method
}//end of class

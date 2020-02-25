package tigge;
import java.util.*;
import java.math.*;
import java.io.*;
import EPSAnalysis.*;

public class Det_MeanError{

        private static final int NUM_TRACKS_MAX = 1000000;

	public static void main(String[] args)throws IOException{
		
		
		String date = null;
		String filename = null;
		
		String fptmp = null;
		
		String[] fpAnalysis = new String[NUM_TRACKS_MAX];
		
		double[] tmpIntDiff = new double[Constants.LT_MAX];
		double[] sumIntDiff = new double[Constants.LT_MAX];
		double[] tmpSepDiff = new double[Constants.LT_MAX];
		double[] sumSepDiff = new double[Constants.LT_MAX];
                double[] tmpIntBias = new double[Constants.LT_MAX];
                double[] sumIntBias = new double[Constants.LT_MAX];
		
		double[] sumIntDiff2 = new double[Constants.LT_MAX];
		double[] sumSepDiff2 = new double[Constants.LT_MAX];
                double[] sumIntBias2 = new double[Constants.LT_MAX];		
		
		
		int[] countDiff = new int[Constants.LT_MAX];
		
		TrackDate forSt;
		
		File ferr=null;
		
		int ntr=0, icur=0;
		
		BufferedWriter writer;
		
		PrintWriter[] det_err = new PrintWriter[NUM_TRACKS_MAX];
		
// Read region file if present

                Region reg = RegionOp.ReadReg();
		
		Arrays.fill(det_err, null);

                if(args.length != 7) {
                  System.err.println("Usage: Det_MeanError [path to dirs] [longitude id] [latitude id] [intensity id] [Match dir] [Mean Out Typ. 0 or 1] [Write Err., 0 or 1]");
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

                System.out.println(forDates.length);
		Arrays.sort(forDates);
		
                FilenameFilter fmatchtyp = new FilenameFilter() {
                   public boolean accept(File d, String name) {
                      boolean ret;
                      if(name != null && (name.startsWith("trmatch") )) {
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
			sumIntDiff[i]=0;
			sumSepDiff[i]=0;
                        sumIntBias[i]=0;
			
			sumIntDiff2[i]=0;
			sumSepDiff2[i]=0;
			sumIntBias2[i]=0;
			
			countDiff[i]=0;
		
		}
		
		for(int i=0; i<forDates.length; i++){
		
			System.out.println(i + " of " + forDates.length);
		        System.out.println(forDates[i].toString());	
			
			date = forDates[i].toString();
                        date = date.substring(date.length() - 10);
				
			System.out.println(date);
				
			forSt = new TrackDate(Integer.parseInt(date.substring(0,4)), Integer.parseInt(date.substring(4,6)), Integer.parseInt(date.substring(6,8)), Integer.parseInt(date.substring(8,10)));

                        File curdir = new File(forDates[i].toString() + "/" + matchd + "/DET/");

                        if (curdir.isDirectory()) {

                            File[] trackFile = curdir.listFiles(fmatchtyp);

                            if(trackFile != null){
                               //System.out.println(trackFile.length);

                               for(int j=0; j<trackFile.length; j++){
                                  for(int l=0; l<65; l++){
                                     tmpSepDiff[l]=-100;
                                     tmpIntDiff[l]=-100;
                                     tmpIntBias[l]=-100;
                                  }
				
		         	  Track analysis = TrackOperations.File2Track(trackFile[j], 0, ilong, ilat, ifd);//analysis track
				  Track forecast = TrackOperations.File2Track(trackFile[j], 1, ilong, ilat, ifd);//forecast track
				  
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
				  
				  
				     if(det_err[icur] == null) {
				        ferr = new File("TRERROR/det_err" + icur + ".dat");
				        det_err[icur] = new PrintWriter(new BufferedWriter(new FileWriter(ferr, true)), true);
				     }
				     
				     det_err[icur].printf("TTT %s %d %n", forSt, 0);
				  
				  }

			  	  tmpIntBias = TrackOperations.intensityDifference2(analysis,forecast,forSt);
                                  tmpIntDiff = TrackOperations.intensityDifference(analysis,forecast,forSt);
				  tmpSepDiff = TrackOperations.separationDistance(analysis,forecast,forSt);
													

				  for(int k=0; k<Constants.LT_MAX; k++){
							
				     if(tmpSepDiff[k] != -100 && tmpIntDiff[k] != -100 && tmpIntBias[k] != -100){
										
				        sumIntDiff[k]+=tmpIntDiff[k];
 					sumIntBias[k]+=tmpIntBias[k];
					sumSepDiff[k]+=tmpSepDiff[k];
					
				        sumIntDiff2[k]+=tmpIntDiff[k] * tmpIntDiff[k];
 					sumIntBias2[k]+=tmpIntBias[k] * tmpIntBias[k];
					sumSepDiff2[k]+=tmpSepDiff[k] * tmpSepDiff[k];
					
					if(iwe == 1) det_err[icur].printf("%4d %e %e %e %n", k, tmpSepDiff[k], tmpIntDiff[k], tmpIntBias[k]);					
								
					countDiff[k]++;
							
				      }

				  }//end of k loop
				  
			       }
			       
			    }
							

			}//end of j for loop
		
		}//end of i for loop


                if(iwe == 1){
                   for(int i=0; i<NUM_TRACKS_MAX; i++){
		      if(det_err[i] != null) det_err[i].close();
		   }
		}

// output deterministic stats

		if(ipt == 0) writer = new BufferedWriter(WriteFile.openfile(null));
		else writer = new BufferedWriter(WriteFile.openfile("det.txt")); 
		
		writer.write("Deterministic Track: Forecast points, Seperation Distance, Intensity Difference, Intensity Bias, Seperation STD, Intensity STD, Intensity Bias STD");
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
					 + Round.round(Math.sqrt(sumIntBias2[i] - sumIntBias[i] * sumIntBias[i]),4)       );
		        writer.newLine();
		        writer.flush();
                }
		
		if(ipt == 1) writer.close();
		
	}//end of main method
}//end of class

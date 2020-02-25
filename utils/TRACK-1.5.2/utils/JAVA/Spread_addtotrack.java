package tigge;
import java.util.*;
import java.io.*;
import EPSAnalysis.*;

public class Spread_addtotrack{
	public static void main(String[] args)throws IOException{
		
		String date = null;
		String filename = null;
		double[] tmpSepDiff = new double[Constants.LT_MAX];
		double[] tmpIntDiff = new double[Constants.LT_MAX];
		double[] sumSepDiff2 = new double[Constants.LT_MAX];
		double[] sumIntDiff2 = new double[Constants.LT_MAX];
		double[] spreadSep = new double[Constants.LT_MAX];
		double[] spreadInt = new double[Constants.LT_MAX];
		
		int[] countDiff2 = new int[Constants.LT_MAX];
		int[] spreadCount= new int[Constants.LT_MAX];
		TrackDate forSt, meanSt, tmpSt;
		File analFile;
		int noOfTracks;
		int ptSt;
		
		boolean flag;
		
	        if(args.length != 5) {
                  System.err.println("Usage: Spread_addtotrack [path to dirs] [longitude id] [latitude id] [intensity id] [Match dir]");
                  System.exit(1);
                }
		
		int ilong = Integer.parseInt(args[1]);
		int ilat = Integer.parseInt(args[2]);
		int ifd = Integer.parseInt(args[3]);

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
				sumSepDiff2[i]=0;
				sumIntDiff2[i]=0;
				countDiff2[i]=0;
				spreadSep[i]=0;
				spreadInt[i]=0;
				spreadCount[i]=0;
				
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
				    //System.out.println(meanfile.toString());
				    
				    if(!meanfile.exists()) continue;

                                    Track mean = TrackOperations.File2Track(meanfile,0, ilong, ilat, ifd);//mean track

                                    for(int l=0; l<mean.getPointNo(); l++){
                                       if(mean.getPoint(l).getIntensity() > 100){
                                          System.out.println(meanfile.toString());
                                       }

                                    }
				    
 // read analysis file

                                    Track analysis = TrackOperations.File2Track(trackFile[j], 0, ilong, ilat, ifd);//analysis track
				    
				    noOfTracks = TrackOperations.getNoOfTracks(trackFile[j]);
				    
				    for(int l=0; l<Constants.LT_MAX; l++){
								
					tmpSepDiff[l]=-100;
					tmpIntDiff[l]=-100;
					sumSepDiff2[l]=0;
					sumIntDiff2[l]=0;
					countDiff2[l]=0;
						
				    }
				    
// open file for writing mean track with individual spreads added

                                    String filestub = trackFile[j].toString().substring(trackFile[j].toString().indexOf("trmatch"));
				    int noOfMeanTracks = TrackOperations.getNoOfTracks(meanfile);
				    if(noOfMeanTracks != 1){
				       System.out.println("Mean track file has too many tracks.");
				       System.exit(1);
				    }

                                    BufferedWriter meanaddspread = new BufferedWriter(new FileWriter(curdir.toString() + "/MEAN/ENSDIFF/" + filestub + "_mean_spread"));
				    
				    //System.out.println(curdir.toString() + "/MEAN/ENSDIFF/" + filestub + "mean_spread");
				    
				    meanSt = mean.getPoint(0).getTrackDate();
				    
				    meanaddspread.write("0\n0 0\nTRACK_NUM       "+ 1 + " ADD_FLD    2   2 &00\n" + "TRACK_ID  " + 0
				    + " START_TIME  " + forSt + "\n");
				    meanaddspread.write("POINT_NUM  " + mean.getPointNo() + "\n");
				    			   
				    ptSt = 0;
                                    tmpSt = new TrackDate(forSt);
				    flag = true;
				    
				    while(flag){
				       if(TrackDate.GreaterThan(meanSt, tmpSt)){
				          ++ ptSt;
					  tmpSt.AddTimeFrames(1);
					  //System.out.println(tmpSt);
				       }
				       else {
				          flag = false;
				       }
	                 	    }
				    
				    
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

					}//end of k loop
						
	
						
				     }//end  of l loop
				     
	
						
				     for(int l=0; l<Constants.LT_MAX; l++){
						
							
					if(countDiff2[l] != 0){
				   	   //System.out.println(countDiff2[l]);
					   spreadSep[l]+=(sumSepDiff2[l]/countDiff2[l]);
					   spreadInt[l]+=(sumIntDiff2[l]/countDiff2[l]);
							
					   spreadCount[l]++;
								
								
								
								
					}
				     }//end of l for loop
				     
				     
				     
				     for(int l=0; l < mean.getPointNo(); l++){
				       
				        if( ! (countDiff2[l + ptSt] == 0)){
				           meanaddspread.write(mean.getPoint(l) + " & " +  String.format("%e", sumSepDiff2[l + ptSt]/countDiff2[l + ptSt]) + " & " + String.format("%e",sumIntDiff2[l + ptSt]/countDiff2[l + ptSt]) + " & \n");
					}
					else {
					   meanaddspread.write(mean.getPoint(l) + " & " +  String.format("%e", 1.0e+25) + " & " + String.format("%e", 1.0e+25) + " & \n");
					} 
				     }
						
				     meanaddspread.close();
	
			       }
			       
			    }		
		
			}//end of j for loop
		
		}//end of i for loop
		

                System.out.println("Spread: Spread Points, Spread Separation, Spread Intensity");
		for(int i=0; i<Constants.LT_MAX; i++){
			System.out.println(Round.round((double)i/4,2) + "\t" + spreadCount[i] + "\t" + Round.round(spreadSep[i]/spreadCount[i],4) + "\t" + Round.round(spreadInt[i]/spreadCount[i],4));
		}
		
		
	}//end of main method
}//end of class

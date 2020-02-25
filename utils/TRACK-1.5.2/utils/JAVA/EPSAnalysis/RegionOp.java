package EPSAnalysis;
import java.util.*;
import java.io.*;

public class RegionOp{

/** Creates a new instance of RegionOp */
    public RegionOp() {
    }


    public static int checkRegion(Region reg, TrackPoint pt){
       int inreg=0;

       double ptlon = pt.getLon(); 
       double ptlat = pt.getLat(); 
       
       if(reg.getLon(0) > reg.getLon(1)){
       
          if(((ptlon >= reg.getLon(0) && ptlon <= 360.0) || (ptlon >= 0.0 && ptlon <= reg.getLon(1))) &&
             (ptlat >= reg.getLat(0) && ptlat <= reg.getLat(1))                                          ) inreg = 1; 
	           
       }
       else {

          if((ptlon >= reg.getLon(0) && ptlon <= reg.getLon(1)) &&
             (ptlat >= reg.getLat(0) && ptlat <= reg.getLat(1))     ) inreg = 1;
	  
       }  
        
       return inreg; 
    } 

    public static Region ReadReg() {
// Read region file if present

       int ireg = 0;

       double lon1, lon2, lat1, lat2;

       Scanner scan=null;

       Region reg=new Region(0.0, 0.0, 0.0, 0.0, 0);

       File file = new File("region.dat");
       try {
          ireg = 1;
          scan = new Scanner(file);

          lon1 = scan.nextDouble();
          lon2 = scan.nextDouble();
          lat1 = scan.nextDouble();
          lat2 = scan.nextDouble();

          reg = new Region(lon1, lon2, lat1, lat2, ireg);

          System.out.println("****INFORMATION****, using region.dat file");

       } catch (FileNotFoundException err) {
          System.out.println("****WARNING****, no region.dat file");
       } finally {
          if (scan != null) scan.close();
       }

       return reg;
    }

}

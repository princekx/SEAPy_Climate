package EPSAnalysis;

public class Region{

    private int ireg;

    private double lon1, lon2;
    private double lat1, lat2; 
    

    Region(double lon1, double lon2, double lat1, double lat2, int ireg){
       this.lon1 = lon1;
       this.lon2 = lon2;
       this.lat1 = lat1;
       this.lat2 = lat2;
       this.ireg = ireg;
    }


    public double getLon(int fl){
	
	return (fl==0) ? lon1 : lon2;
	
    }	
	
    public double getLat(int fl){
	
	return (fl==0) ? lat1 : lat2;
	
    }

    public int getIreg(){

        return ireg;

    }

}

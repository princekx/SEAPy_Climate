package EPSAnalysis;
import java.io.*;

public class WriteFile{


    public static Writer openfile(String fileName) throws IOException, FileNotFoundException {

           if (fileName != null) {
               return new PrintWriter(fileName);
           }
           else
               return new OutputStreamWriter(System.out);
    }



}

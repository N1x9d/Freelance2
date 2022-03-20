
import java.util.concurrent.CountDownLatch;

//import java.util.Arrays;

public class NewThreadH implements Runnable {

    float[][] Ez;
    float[][] Hx;
    float[][] Hy;
    int iBegin;
    int iEnd;
    float tt;
    int nt = 1;
    
    CountDownLatch cdl;

    public NewThreadH(float[][] E, float[][] H, float[][] H2,
            int iBegi, int iEn, float ttt, CountDownLatch cdl) {
        //Ez = Arrays.copyOf(E, E.length);
        //Hx = Arrays.copyOf(H, H.length);
        //Hy = Arrays.copyOf(H2, H2.length);
        this.cdl = cdl;
        
        Ez = E;
        Hx = H;
        Hy = H2;
        iBegin = iBegi;
        iEnd = iEn;
        tt = ttt;
      
    }

    @Override
    public void run() {
        
        //__ магнитное поле перенормирвано для облегчения рассчета
        for (int i = iBegin; i <= iEnd; i++) {          
            for (int j = 1; j <= tt; j++) {
                Hx[i][j] += Ez[i][j] - Ez[i][j + 1];
                Hy[i][j] += Ez[i + 1][j] - Ez[i][j];
            }
        }

        cdl.countDown();
    }
}

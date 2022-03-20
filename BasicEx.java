
import com.meapsoft.FFT;
import java.awt.Color;
import java.awt.EventQueue;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import static java.lang.Math.*;
import java.util.Arrays;
import java.util.Calendar;
import javax.swing.JFrame;
import javax.swing.JPanel;



public class BasicEx extends JFrame {

    public static int nx = 2048;//limit=15M/1gb одной компоненты. ПОВЫСИЛ лимит памяти, и уже Флоат, но 2 компоненты
    public static int ny = 500;//2048*2000(y) = 104sec 2компоненты
    public static double[][] forFurier = new double[2][1];

    public static float[][]Hx;
    public static float[][]Hy;
    public static float[][]Ez;

    private static String[] arg;

    public BasicEx() {
        if (arg.length > 1) {
            nx = Integer.parseInt(arg[0]);
            ny = Integer.parseInt(arg[1]);
        }
    // начало замера
        long tTime = Calendar.getInstance().getTimeInMillis();
        NewThread co = new NewThread(nx, ny, "cos");
        NewThread si = new NewThread(nx, ny, "sin");
        try {
            co.t.join();
            si.t.join();
        } catch (InterruptedException e) {
        }
        System.out.println("Full calculate time   "
                + (double) (Calendar.getInstance().getTimeInMillis() - tTime) / 1000 + " sec");
        FDTD.service.shutdown();
        //
        for (int i = 0; i < nx + 1; i++) {
            for (int j = 0; j < ny + 1; j++) {
                co.E[i][j] = co.E[i][j] * co.E[i][j] + si.E[i][j] * si.E[i][j];
            }
        }
        si = null;
        // ТЕПЕРЬ ЭТО ИНТЕНСИВНОсТЬ

        //FFT fft = new FFT(nx);
        //fft.fft(forFurier[0], forFurier[1]);
        //for (int i = 0; i < forFurier[0].length; i++) {
//            forFurier[0][i] = pow( //amplitude //in low power
  //                  pow(forFurier[0][i], 2) + pow(forFurier[1][i], 2), 0.5);

    //    }
        //System.out.println(Arrays.toString(forFurier[0]));
        ny = 1366;
        //add(new Surface(co.E));
        //setTitle("FDTD Java 2D example");
        //setSize(ny, nx + 40);
        //setLocationRelativeTo(null);
        //setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    }

    public static void main(String[] args) {
        arg = args;

        EventQueue.invokeLater(new Runnable() {

            @Override
            public void run() {
                BasicEx ex = new BasicEx();
                ex.setVisible(true);
            }
        });
    }

   /* class Surface extends JPanel {

    float[][] arr = new float[1][1];

    public Surface(float[][] arr) {
        super();
        this.arr = arr; // uge ABS

    }

    @Override
    public void paintComponent(Graphics g) {
        super.paintComponent(g);

        int x = arr.length;
        int y = arr[0].length;
        double mAx = 0;
        for (int i = 0; i < x; i++) {
            for (int j = 0; j < y; j++) {
                mAx = max(mAx, arr[i][j]);
            }
        }

        BufferedImage bufferedImage = new BufferedImage(y, x, BufferedImage.TYPE_INT_RGB);
        int rgb;
        for (int i = 0; i < x; i++) {
            for (int j = 0; j < y; j++) {
                rgb = (int) (arr[i][j] / mAx * 255 * 3);
                if (rgb < 255) {
                    rgb = (new Color(rgb, 0, 0)).hashCode();
                } else if (rgb < 510) {
                    rgb = (new Color(255, (rgb - 255), 0)).hashCode();
                } else {
                    rgb = (new Color(255, 255, (rgb - 510))).hashCode(); //rgb =  (new Color(rgb, rgb, rgb)).hashCode();
                }
                bufferedImage.setRGB(j, i, rgb);//(x=j, y=i)
            }
        }
        Graphics2D g2d = (Graphics2D) g;
        g2d.drawImage(bufferedImage, null, this);

        double[] f = BasicEx.forFurier[0];
        int first = max(BasicEx.ny - BasicEx.nx, 0);
        mAx = Arrays.stream(f).max().getAsDouble();
        g2d.setColor(Color.red);
        g2d.drawString("maximal spectr =  " + Double.toString(mAx), 1100, 100);
        x = BasicEx.nx / 2;
        int tX = min(700, BasicEx.nx);
        int shift;
        int i;
        //Antialiasing ON
        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
                RenderingHints.VALUE_ANTIALIAS_ON);
        for (i = 0; i < BasicEx.nx - 1; i++) {
            shift = (i < x) ? x : -x;
            g2d.drawLine(first + i + shift, (int) ((mAx - f[i]) / mAx * tX),
                    first + i + 1 + shift, (int) ((mAx - f[i + 1]) / mAx * tX));
        }
        g2d.drawLine(first + x - 1, (int) ((mAx - f[i]) / mAx * tX),
                first + x, (int) ((mAx - f[0]) / mAx * tX));
    }*/
//}
}

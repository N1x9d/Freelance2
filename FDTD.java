import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import static java.lang.Math.*;


public class FDTD {

    public static ExecutorService service = Executors.newFixedThreadPool(4);// fixed threads

    public static float[][] callFDTD(int nx, int ny, String method) {
        CountDownLatch cdl;

        int i, j;//for loops
        float x; //for coordinate
        final float lambd = 1064e-9f; // Длина волны в метрах
        final float dx = lambd / 15;         // размер шага сетки. Максимум λ/10. Чем меньше, тем выше точность.
        final float dy = dx; //реально в коде не используется. Ячейки квадратные
        final float period = 2e-6f; // поперечный период решетки фотонного кристалла
        final float Q = 1.08f;//1.0 = default // специфическое отношение
        final float n = 1;// минимальный показатель преломления кристалла

        final float prodol = 2 * n * period * period / lambd / Q; //продольный период
        // основные массивы полей
        float[][] Ez = new float[nx + 1][ny + 1];
        float[][] Hx = new float[nx - 1 + 1][ny - 1 + 1];
        float[][] Hy = new float[nx][ny];

        final float omega = (float) (2 * PI / lambd); // циклическая частота световой волны
        final float dt = dx / 2; // шаг по времени. НЕ больше, чем по координате
        final float tau = 2e-5f * 999;// время затухания источника света. Нужно для запуска коротких импульсов. Если большое, затухание незаметно.
        final float s = dt / dx; //отношение шага по времени и пространству
        final float k3 = (1 - s) / (1 + s);// для граничных условий
        final float w = 19e-7f;// полуширина гауссова пучка
        final float alpha = (float) (sin(0.0 / 180 * PI));//radians. // угол лазерного пучка от нормали. В радианах.
        final float Z = 376.7303f;



        final int begin = 10;// координата, где пустота сменяется кристаллом. Без необходимости много не ставить.
        System.out.println("Imeem__Sloy= " + (ny - begin) * dy / prodol * 2);
        final float mod = 0.008f * 2;//максимальная глубина модуляции диэлектрической проницаемости = 2*Δn;

        final float ds = dt * dt / dx / dx;// renormalization  dH=dt/dx/Z; // константа для подгона систем исчисления
        float[][] e = new float[nx + 1][ny + 1];
        //Инициализация диэлектрической проницаемости (решетка кристалла)
        for (i = 1; i < nx + 1; i++) {
            for (j = 1; j < ny + 1; j++) {
                // EPSILON MATRIX full initialization
                e[i][j] = (float) (ds / (n + ((j < begin) ? 0 : (mod / 2) * (1 + signum(-0.1 
                           + cos(2 * PI * (i - nx / 2.0 + 0.5) * dx / period) * sin(2 * PI * (j - begin) * dy / prodol))))));
            }
        }
        //Массивы, используемые для граничных условий:
        float[][] end = new float[2][nx + 1]; // boundary conditions
        float[][] top = new float[2][ny + 1];
        float[][] bottom = new float[2][ny + 1];
        //Предел счета по времени
        final int tmax = (int) (ny * 2.2);
        System.out.println("START CICLE");
        for (int t = 1; t <= tmax; t++) {                   
            float tt = Math.min(t * s + 10, ny - 1);
            //gauss 
            switch (method) {
                case "cos":        
                    for (i = 1; i <= nx - 1; i++) {
                        x = (float) (dx * (i - (float) nx / 2 + 0.5)); // поперечная координата в метрах
                        Ez[i][1] = (float) (exp(-pow(x, 2) / w / w - (t - 1) * dt / tau) * cos((x * alpha + (t - 1) * dt) * omega));//перепишешь лямбді массива
                    }
                    break;
                case "sin":
                    for (i = 1; i <= nx - 1; i++) {
                        x = (float) (dx * (i - (float) nx / 2 + 0.5)); // поперечная координата в метрах
                        Ez[i][1] = (float) (exp(-pow(x, 2) / w / w - (t - 1) * dt / tau) * sin((x * alpha + (t - 1) * dt) * omega));
                    }
                    break;
            }

            for (i = 1; i <= nx; i++) {  // boundary conditions
                end[0][i] = Ez[i][ny - 1];
                end[1][i] = Ez[i][ny];
            }
            //Дальше копируем временные массивы под поглощающие граничные условия Мура:
            System.arraycopy(Ez[1], 0, top[0], 0, ny + 1);
            System.arraycopy(Ez[2], 0, top[1], 0, ny + 1);
            System.arraycopy(Ez[nx - 1], 0, bottom[0], 0, ny + 1);
            System.arraycopy(Ez[nx], 0, bottom[1], 0, ny + 1);

            // main Ez
            // расчет EZ в параллельных потоках
            cdl = new CountDownLatch(2);
            NewThreadE first = new NewThreadE(Ez, Hx, Hy, e, 2, nx / 2, tt, cdl);
            NewThreadE second = new NewThreadE(Ez, Hx, Hy, e, nx / 2 + 1, nx - 1, tt, cdl);
            service.execute(first);
            service.execute(second);
            try {
                cdl.await();
            } catch (InterruptedException ex) {
            }

            for (i = 1; i <= nx; i++) {    // boundary conditions
                Ez[i][ny] = end[0][i] + k3 * (end[1][i] - Ez[i][ny - 1]);//end
            }
            for (i = 1; i <= ny; i++) {
                Ez[1][i] = top[1][i] + k3 * (top[0][i] - Ez[2][i]);//verh kray
                Ez[nx][i] = bottom[0][i] + k3 * (bottom[1][i] - Ez[nx - 1][i]);
            }
            //Лазер
            switch (method) {
                case "cos":
                    for (i = 1; i <= nx - 1; i++) {
                        x = (float) (dx * (i - (float) nx / 2 + 0.5));
                        Ez[i][1] = (float) (exp(-pow(x, 2) / w / w - (t - 1) * dt / tau) * cos((x * alpha + t * dt) * omega));
                    }
                    break;
                case "sin":
                    for (i = 1; i <= nx - 1; i++) {
                        x = (float) (dx * (i - (float) nx / 2 + 0.5));
                        Ez[i][1] = (float) (exp(-pow(x, 2) / w / w - (t - 1) * dt / tau) * sin((x * alpha + t * dt) * omega));
                    }
                    break;
            }

            cdl = new CountDownLatch(2);
            //магнитное поле многопоток
            NewThreadH firstH = new NewThreadH(Ez, Hx, Hy, 1, nx / 2, tt, cdl);
            NewThreadH secondH = new NewThreadH(Ez, Hx, Hy, nx/2+1, nx-1, tt, cdl);
            service.execute(firstH);
            service.execute(secondH);
            try {
                cdl.await();
            } catch (InterruptedException ex) {
            }

        }

        int pos = method.equals("cos") ? 0 : 1;//какую из компонент записывать
        //BasicEx.forFurier[pos] = new double[nx];
        int endF = (int) (ny * 0.95);//0.99
        for (i = 1; i <= nx; i++) {
            //BasicEx.forFurier[pos][i - 1] = Ez[i][endF];
            for (j = 1; j <= ny; j++) {
                Ez[i][j] = abs(Ez[i][j]);// amplitude
            }
        }
        BasicEx.Hx=Hx;
        BasicEx.Hy=Hy;
        BasicEx.Ez=Ez;
        Hx = null;
        Hy = null;
        e = null;
        return Ez;

    }

}


public class NewThread implements  Runnable {
    Thread  t; 
    private int nx=0, ny=0;
    String met = "cos";
    float [][]E;
    
    public NewThread(int x, int y, String method) {
        nx = x;
        ny = y;
        met = method;
        t = new Thread (this, "Next Thread");
        t.start();
        t.setPriority(10);
    }
    
    @Override
    public  void  run ()  {
            E =  FDTD.callFDTD(nx, ny, met);;
        }
}

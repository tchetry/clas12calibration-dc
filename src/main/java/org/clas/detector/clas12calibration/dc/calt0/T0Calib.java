/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clas.detector.clas12calibration.dc.calt0;
import org.clas.detector.clas12calibration.dc.t2d.TableLoader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.clas.detector.clas12calibration.dc.analysis.Coordinate;
import org.clas.detector.clas12calibration.dc.calt2d.SegmentProperty;
import org.clas.detector.clas12calibration.viewer.AnalysisMonitor;
import org.jlab.detector.calib.utils.CalibrationConstants;
import org.jlab.groot.data.H1F;
import org.jlab.groot.group.DataGroup;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent; 
import org.jlab.detector.calib.utils.ConstantsManager;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.math.F1D;
import org.jlab.io.hipo.HipoDataSource;
import org.jlab.jnp.hipo4.data.SchemaFactory;
import org.jlab.rec.dc.Constants;
import org.jlab.rec.dc.hit.FittedHit;
import org.jlab.utils.groups.IndexedList;
import org.jlab.utils.system.ClasUtilsFile;
/**
 *
 * @author KPAdhikari, ziegler
 */
public class T0Calib extends AnalysisMonitor{
    //public HipoDataSync writer = null;
    //private HipoDataEvent hipoEvent = null;
    private SchemaFactory schemaFactory = new SchemaFactory();
    PrintWriter pw = null;
    File outfile = null;
    private int runNumber;
    private String analTabs = "Corrected TDC";;
    public T0Calib(String name, ConstantsManager ccdb) throws FileNotFoundException {
        super(name, ccdb);
        this.setAnalysisTabNames(analTabs);
        this.init(false, "T0");
        outfile = new File("Files/ccdbConstantst0.txt");
        pw = new PrintWriter(outfile);
        pw.printf("#& Sector Superlayer Slot Cable T0Correction T0Error\n");
        
        String dir = ClasUtilsFile.getResourceDir("CLAS12DIR", "etc/bankdefs/hipo4");
        schemaFactory.initFromDirectory(dir);
       
        if(schemaFactory.hasSchema("TimeBasedTrkg::TBHits")) {
            System.out.println(" BANK FOUND........");
        } else {
            System.out.println(" BANK NOT FOUND........");
        }
        //writer = new HipoDataSync(schemaFactory);
        //writer.setCompressionType(2);
        //hipoEvent = (HipoDataEvent) writer.createEvent();
        //writer.open("TestOutPut.hipo");
        //writer.writeEvent(hipoEvent);
        
        
        
    }
    int nsl  = 6;
    int nsec = 6;
    int nCrates = 18;// Goes from 41 to 58 (one per chamber)
    int nSlots = 20; // Total slots in each crate (only 14 used)
    int nChannels = 96;// Total channels per Slot (one channel per wire)
    int nLayers0to35 = 36;// Layers in each sector (0th is closest to CLAS
    int nCables = 84;
    int nCables6 = 6; // # of Cables per DCRB or STB.
    int nSlots7 = 7; // # of STBs or occupied DCRB slots per SL.
    boolean[][][][] Fitted = new boolean[nsec][nsl][nSlots7][nCables6];
    int[] nTdcBins =
    { 50, 50, 50, 50, 50, 50 };
    int[] nTimeBins =
    { 50, 50, 50, 50, 50, 50 };
    double[] tLow =
    { 80.0, 80.0, 80.0, 80.0, 80.0, 80.0 };
    
    public static final double[] tLow4T0Fits  = {-40.0, -40.0, -40.0, -40.0, -40.0, -40.0};
    public static final double[] tHigh4T0Fits  = {380.0, 380.0, 680.0, 780.0, 1080.0, 1080.0}; 

    public static  double[][][][] fitMax ;


    //H1F[][][][] h = new H1F[6][6][nSlots7][nCables6];
    private Map<Coordinate, H1F> TDCHis        = new HashMap<Coordinate, H1F>();    
    public  Map<Coordinate, FitLine> TDCFits   = new HashMap<Coordinate, FitLine>();
    public  Map<Coordinate, Double> T0s        = new HashMap<Coordinate, Double>();
    
    @Override
    public void createHistos() {
        //histo max range for the fit
        fitMax = new double[nsec][nsl][nSlots][nCables]; 
        // initialize canvas and create histograms
        this.setNumberOfEvents(0);
        DataGroup hgrps = new DataGroup(6,7);
        String hNm;
        String hTtl;
        int ijk=-1;
        for (int i = 0; i < nsec; i++)
        {
            for (int j = 0; j < nsl; j++)
            {
                for (int k = 0; k < nSlots7; k++)
                {
                    for (int l = 0; l < nCables6; l++)
                    {
                        hNm = String.format("timeS%dS%dS%dCbl%d", i + 1, j + 1, k + 1, l + 1);
                        
                        TDCHis.put(new Coordinate(i,j,k, l), new H1F(hNm, 150, tLow4T0Fits[j], tHigh4T0Fits[j])); 
                                                                                                                                                                                        // HBHits
                        hTtl = String.format("time (Sec%d SL%d Slot%d Cable%d)", i + 1, j + 1, k + 1, l + 1);
                        TDCHis.get(new Coordinate(i,j,k, l)).setTitleX(hTtl);
                        TDCHis.get(new Coordinate(i,j,k, l)).setLineColor(1);
                        TDCFits.put(new Coordinate(i,j,k, l), new FitLine());
                        hgrps.addDataSet(TDCHis.get(new Coordinate(i, j, k, l)), 0);
                        
                        T0s.put(new Coordinate(i,j,k, l), 0.0);
                        
                        Fitted[i][j][k][l] = false;
                    }
                    this.getDataGroup().add(hgrps, i+1, j+1, k+1);
                }
                
            }
        }

        this.getDataGroup().add(hgrps,0,0,0);
        for (int i = 0; i < nsec; i++) {
            for (int j = 0; j < nsl; j++) {
                for (int k = 0; k < nSlots7; k++) {
                    this.getCalib().addEntry(i+1,j+1,k+1);
                }
            }
        }
        this.getCalib().fireTableDataChanged();
    }
     
    @Override
    public void plotHistos() {
        this.getAnalysisCanvas().getCanvas(analTabs).setGridX(false);
        this.getAnalysisCanvas().getCanvas(analTabs).setGridY(false);
        this.getAnalysisCanvas().getCanvas(analTabs).divide(nCables6/2, 2);
        this.getAnalysisCanvas().getCanvas(analTabs).update();
        
        
    }
    @Override
    public void timerUpdate() {
    }
    
    @Override
    public void analysis() {
        this.plotFits();
    }
    public void plotFits() {
        File file2 = new File("");
	file2 = outfile;
        DateFormat df = new SimpleDateFormat("MM-dd-yyyy_hh.mm.ss_aa");
        String fileName = "Files/ccdb_run" + this.runNumber + "time_" 
                + df.format(new Date())  + ".txt";
        file2.renameTo(new File(fileName));
        
        for (int i = 0; i < nsec; i++)
        {
            for (int j = 0; j < nsl; j++)
            {
                for (int k = 0; k < nSlots7; k++)
                {
                    for (int l = 0; l < nCables6; l++)
                    { 
                        // if(this.fitThisHisto(this.TDCHis.get(new Coordinate(i,j,k,l)))==true) {
                            this.runFit(i, j, k, l);
                            int binmax = this.TDCHis.get(new Coordinate(i,j,k,l)).getMaximumBin();
                            fitMax[i][j][k][l] = this.TDCHis.get(new Coordinate(i,j,k,l)).getDataX(binmax);
                            
                        // }
                    }
                }
            }
        }
        pw.close();
        this.getCalib().fireTableDataChanged();  
        
    }
    
    public int NbRunFit = 0;
    public void runFit(int i, int j, int k, int l) {
        
        System.out.println(" **************** ");
        System.out.println(" RUNNING THE FITS ");
        System.out.println(" **************** "); 
	
        double[] Tminmax = this.getT0(i, j, k, l);
        
        //Sector Superlayer Slot Cable T0Correction T0Error
        pw.printf("%d\t %d\t %d\t %d\t %.6f\t %.6f\n",
            (i+1), (j+1), (k+1), (l+1), 
            Tminmax[0], 
            Tminmax[1]);
        System.out.printf("%d\t %d\t %d\t %d\t %.6f\t %.6f\n",
            (i+1), (j+1), (k+1), (l+1), 
            Tminmax[0], 
            Tminmax[1]);
        
        Fitted[i][j][k][l] = true;
        System.out.println(" FITTED ? "+Fitted[i][j][k][l]);
    }
     private void updateTable(int i, int j,  int k, double t0) {
       this.getCalib().setDoubleValue(t0, "T0", i+1, j+1, k+1);
     }
    
    int counter = 0;
    public  HipoDataSource reader = new HipoDataSource();
    

    int count = 0;
    public static int polarity =-1;
    public List<FittedHit> hits = new ArrayList<>();
    Map<Integer, ArrayList<Integer>> segMapTBHits = new HashMap<Integer, ArrayList<Integer>>();
    Map<Integer, SegmentProperty> segPropMap = new HashMap<Integer, SegmentProperty>();
    //List<FittedHit> hitlist = new ArrayList<>();
    private ReadTT cableMap = new ReadTT();
    @Override
    public void processEvent(DataEvent event) {
        
        if (!event.hasBank("RUN::config")) {
            return ;
        }
        
        DataBank bank = event.getBank("RUN::config");
        int newRun = bank.getInt("run", 0);
        if (newRun == 0) {
           return ;
        } else {
           count++;
        }
        
        if(count==1) {
            Constants.Load();
            TableLoader.FillT0Tables(newRun, "default");
            ReadTT.Load(newRun, "default"); 
            runNumber = newRun; 
        }
        if(!event.hasBank("TimeBasedTrkg::TBHits")) {
            return;
        } 
        // get segment property
        
        DataBank bnkHits = event.getBank("TimeBasedTrkg::TBHits");
        
        for (int j = 0; j < bnkHits.rows(); j++) {
            
            int sec = bnkHits.getInt("sector", j);
            int sl = bnkHits.getInt("superlayer", j);
            int lay = bnkHits.getInt("layer", j);// layer goes from 1 to 6 in data
            int wire = bnkHits.getInt("wire", j);// wire goes from 1 to 112 in data
            //int lay0to35 = (sl - 1) * 6 + lay - 1;
            //int region0to2 = (int) ((sl - 1) / 2);
            int slot1to7  = (int) ((wire - 1) / 16) + 1;
            int wire1to16 = (int) ((wire - 1) % 16 + 1);
            int cable1to6 = ReadTT.CableID[lay - 1][wire1to16 - 1];
            double time = (double) bnkHits.getFloat("time", j)
                    + (double) bnkHits.getFloat("T0", j);   
            if(bnkHits.getByte("trkID", j)!=-1)
                this.TDCHis.get(new Coordinate(sec-1, sl-1, slot1to7-1, cable1to6-1))
                    .fill(time);
        } 
    }
    
    public void Plot(int i , int j, int k) {
        
        for (int l = 0; l < nCables6; l++){
            this.getAnalysisCanvas().getCanvas(analTabs).cd(l);
            this.getAnalysisCanvas().getCanvas(analTabs)
                    .draw(this.TDCHis.get(new Coordinate(i, j, k, l)));
            
            if(Fitted[i][j][k][l]==true) {
                this.getAnalysisCanvas().getCanvas(analTabs).cd(l);
                            this.getAnalysisCanvas().getCanvas(analTabs)
                                .draw(this.TDCFits.get(new Coordinate(i, j, k, l)), "same");
            }
        }
    }
    
    @Override
    public void constantsEvent(CalibrationConstants cc, int col, int row) {
        String str_sector    = (String) cc.getValueAt(row, 0);
        String str_layer     = (String) cc.getValueAt(row, 1);
        String str_slot     = (String) cc.getValueAt(row, 2);
        System.out.println(str_sector + " " + str_layer + " " );
        IndexedList<DataGroup> group = this.getDataGroup();

       int sector    = Integer.parseInt(str_sector);
       int layer     = Integer.parseInt(str_layer);
       int slot      = Integer.parseInt(str_slot);

       if(group.hasItem(sector,layer,slot)==true){
           this.Plot(sector-1, layer-1, slot-1);
       } else {
           System.out.println(" ERROR: can not find the data group");
       }
   
    }

    private boolean fitThisHisto(H1F h) {
        boolean pass = false;
        int nevent = 0;
        int maxbin = h.getMaximumBin();
        for (int ix =0; ix< maxbin; ix++) {
            double y = h.getBinContent(ix);
            double err = h.getBinError(ix);
            
            if(err>0 && y>0) {
                nevent+=y;
            }
        }
        
        if(nevent > 99)
            pass = true;
        
        return pass;
    }

    private double getThreshold(H1F h) {
        // find the bin at which the integral corresponds to 1% of the full integral to max 
        // this is the max range to obtain a threshold from a flat line fit
        double integral = 0;
        double partintegral = 0;
        
        GraphErrors gr = new GraphErrors(); 
        for (int ix =0; ix< h.getMaximumBin(); ix++) {
            integral+= h.getBinContent(ix);
        }
        double x = 0;
        for (int ix =0; ix< h.getMaximumBin(); ix++) {
            x = h.getDataX(ix);
            double y = h.getBinContent(ix);
            double err = h.getBinError(ix);
            if(err<1) {
                err = 1.4142;
            }
            // fill graph
            gr.addPoint(x, y, 0, err);
            partintegral += h.getBinContent(ix);
            
            if(partintegral>0.5*integral)
                break;
        }
        // fit the graph 
        F1D f0 = new F1D("f0","[p0]", h.getDataX(0), x);
        f0.setParameter(0, 0);
        DataFitter.fit(f0, gr, "Q"); 
        return f0.getParameter(0);
    }
    private double[] getT0(int i, int j, int k, int l) {
        System.out.println("Getting t0 for i,j,k,l = "+i+" "+j+" "+k+" "+l );
        H1F h = this.TDCHis.get(new Coordinate(i,j,k,l));
        
        double thres = 0;//this.getThreshold(h);
        double [] T0val = new double[2];
        F1D f1 = new F1D("f1","[a]*x+[b]", h.getDataX(0), h.getDataX(20));
        
        F1D gausFunc = new F1D("gausFunc", "[amp]*gaus(x,[mean],[sigma])+[p0]", 
                h.getDataX(0), h.getDataX(h.getMaximumBin())); 
        
        gausFunc.setParameter(0, h.getMax());
        gausFunc.setParameter(1, -0.0);
        gausFunc.setParameter(2, 0.05);
        gausFunc.setParameter(3, 0);
        
        DataFitter.fit(gausFunc, h, "Q"); 
        
        double tmidY = gausFunc.getParameter(0)/2;
        double tminY = gausFunc.getParameter(3);
        double del_min_halfmaxY = tmidY-tminY;
        
        double minRangeY = tmidY-del_min_halfmaxY/2;
        double maxRangeY = tmidY;
        if(h.getMax()>tmidY && tmidY+(h.getMax()-tmidY)/3 < h.getMax() ) {
            maxRangeY+=(h.getMax()-tmidY)/3;
        }
        
        System.out.println(" minRangeY "+minRangeY+" maxRangeY "+maxRangeY);
        
        f1.setParameter(0, 0);
        f1.setParameter(1, 0);

        GraphErrors gr = new GraphErrors(); 
        
        int t0idx  = -1;
        int t0midx = -1;
        double t0 = Double.NEGATIVE_INFINITY; 
        for (int ix =0; ix< h.getMaximumBin(); ix++) {
            if(h.getBinContent(ix)>=maxRangeY) {
                t0midx= ix;
                break;
            }
        }
        for (int ix =0; ix< h.getMaximumBin(); ix++) {
            if(h.getBinContent(ix)>=minRangeY) {
                t0idx= ix;
                break;
            }
        }
        int diffBins = t0midx - t0idx;
        System.out.println("diffBins "+diffBins);
//        for (int ix =0; ix< h.getMaximumBin(); ix++) {
//            if(h.getBinContent(ix) >thres 
//                        && t0 == Double.NEGATIVE_INFINITY) {
//                    t0 = h.getDataX(ix);
//                    t0idx = ix;
//                break;
//            }
//        }
        for (int ix =t0idx; ix< t0midx; ix++) {
            gr.addPoint(h.getDataX(ix), h.getBinContent(ix), 0, h.getBinError(ix));
        }
        
        if(gr.getDataSize(0)>1) {
            f1.setRange(h.getDataX(t0idx), h.getDataX(t0midx));
            DataFitter.fit(f1, gr, "Q"); 
        }
        
        double n = tminY-f1.getParameter(1);
        double d = f1.getParameter(0);
        double en = -f1.parameter(1).error();
        double ed = f1.parameter(0).error();
        double T0 = n/d;
        double T0Err = this.calcError(n, en, d, ed);
	if(Double.isNaN(T0)|| Double.isNaN(T0Err)){
            T0 = 0.0;
            T0Err = 0.0;
        }
        T0val[1] =T0Err;
        T0val[0] = T0;
        h.setOptStat(0);
        String t = "T0 = "+(float)T0;
        h.setTitle(t);
        T0s.put(new Coordinate(i,j,k, l), T0);
        this.updateTable(i, j, k, T0);
        TDCFits.put(new Coordinate(i,j,k,l), 
                new FitLine("f"+""+i+""+j+""+k+""+l, i, j, k, l, 
                T0, h.getDataX(t0midx+diffBins/2)) );
        TDCFits.get(new Coordinate(i,j,k,l)).setLineStyle(4);
        TDCFits.get(new Coordinate(i,j,k,l)).setLineWidth(5);
        TDCFits.get(new Coordinate(i,j,k,l)).setLineColor(8);
        TDCFits.get(new Coordinate(i,j,k,l)).setParameters(new double[] {f1.getParameter(0), f1.getParameter(1)});
        
        return T0val;
    }

    private double calcError(double n, double en, double d, double ed) {
        return Math.sqrt((en/d)*(en/d)+(ed*n/(d*d))*(ed*n/(d*d)));
    }
}


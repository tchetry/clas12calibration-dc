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
import org.clas.detector.clas12calibration.dc.calt0.FitFunction;
import org.clas.detector.clas12calibration.dc.calt0.FitLine;
import org.clas.detector.clas12calibration.dc.calt2d.SegmentProperty;
import org.clas.detector.clas12calibration.viewer.AnalysisMonitor;
import org.freehep.math.minuit.FCNBase;
import org.freehep.math.minuit.FunctionMinimum;
import org.freehep.math.minuit.MnMigrad;
import org.freehep.math.minuit.MnScan;
import org.freehep.math.minuit.MnUserParameters;
import org.jlab.detector.calib.utils.CalibrationConstants;
import org.jlab.groot.data.H1F;
import org.jlab.groot.group.DataGroup;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent; 
import org.jlab.detector.calib.utils.ConstantsManager;
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
    private int runNumber;
    private String analTabs = "Corrected TDC";;
    public T0Calib(String name, ConstantsManager ccdb) throws FileNotFoundException {
        super(name, ccdb);
         this.setAnalysisTabNames(analTabs);
        this.init(false, "T0");
        File outfile = new File("Files/ccdbConstantst0.txt");
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
    private Map<Coordinate, H1F> TDCHis                     = new HashMap<Coordinate, H1F>();    
    private Map<Coordinate, FitFunction> TDCFit             = new HashMap<Coordinate, FitFunction>();
    private Map<Coordinate, MnUserParameters> TDCFitPars    = new HashMap<Coordinate, MnUserParameters>();
    public  Map<Coordinate, FitLine> TDCFits                = new HashMap<Coordinate, FitLine>();
    
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
                        TDCFits.put(new Coordinate(i,j,k, l), new FitLine());;
                        hgrps.addDataSet(TDCHis.get(new Coordinate(i, j, k, l)), 0);
                        
                        
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
        
        pw.close();
        File file2 = new File("");

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
                        this.runFit(i, j, k, l);
                        int binmax = this.TDCHis.get(new Coordinate(i,j,k,l)).getMaximumBin();
                        fitMax[i][j][k][l] = this.TDCHis.get(new Coordinate(i,j,k,l)).getDataX(binmax);
                        TDCFits.put(new Coordinate(i,j,k,l), new FitLine("f"+""+i+""+j+""+k+""+l, i, j, k, l, tLow4T0Fits[j], fitMax[i][j][k][l] , 
                        TDCFitPars.get(new Coordinate(i,j,k,l))));
                        TDCFits.get(new Coordinate(i,j,k,l)).setLineStyle(4);
                        TDCFits.get(new Coordinate(i,j,k,l)).setLineWidth(5);
                        TDCFits.get(new Coordinate(i,j,k,l)).setLineColor(8);
                        
                    }
                }
            }
        }
            
        this.getCalib().fireTableDataChanged();  
        
    }
    private int maxIter = 10;
    
    private MnScan  scanner = null;
    private MnMigrad migrad = null;
    
    public int NbRunFit = 0;
    public void runFit(int i, int j, int k, int l) {
        
        System.out.println(" **************** ");
        System.out.println(" RUNNING THE FITS ");
        System.out.println(" **************** ");
        TDCFit.put(new Coordinate(i, j, k, l), 
                new FitFunction(TDCHis.get(new Coordinate(i, j, k, l))));
        
        scanner = new MnScan((FCNBase) TDCFit.get(new Coordinate(i, j, k, l)), 
                TDCFitPars.get(new Coordinate(i, j, k, l)),2);
	//scanner.fix(2);
        System.out.println(" MINIMIZING............. ");
        FunctionMinimum scanmin = scanner.minimize();
        if(scanmin.isValid())
            TDCFitPars.put(new Coordinate(i, j, k, l),scanmin.userParameters());
        
        migrad = new MnMigrad((FCNBase) TDCFit.get(new Coordinate(i, j, k, l)), 
                TDCFitPars.get(new Coordinate(i, j, k, l)),1);
        migrad.setCheckAnalyticalDerivatives(true);
        
        FunctionMinimum min = null ;
        
        
        for(int it = 0; it<maxIter; it++) {
            
            min = migrad.minimize();
            System.err.println("****************************************************");
            System.err.println("*   FIT RESULTS  FOR SUPERLAYER  "+(j+1)+" at iteration "+(it+1)+"  *");
            System.err.println("****************************************************");  
            
            if(min.isValid()) {
                TDCFitPars.put(new Coordinate(i, j, k, l),min.userParameters());  
            }
            
            if(min!=null)
                System.err.println(min);
        }
            
        //Sector Superlayer Slot Cable T0Correction T0Error
        pw.printf("%d\t %d\t %d\t %d\t %.6f\t %.6f\t "
                + "%.6f\t %.6f\t %.6f\t %.6f\t %.6f\t %d\t %.6f\t %.6f\t %d\n",
            (i+1), (j+1), (k+1), (l+1), 
            TDCFitPars.get(new Coordinate(i, j, k, l)).value(0), 
            TDCFitPars.get(new Coordinate(i, j, k, l)).error(0));
        
        
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
            this.loadFitPars(); 
            polarity = (int)Math.signum(event.getBank("RUN::config").getFloat("torus",0));
            runNumber = newRun;
            System.out.println("CONSTANTS LOADED!!!!!!!!!!!!");
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
            
            this.TDCHis.get(new Coordinate(sec-1, sl-1, slot1to7-1, cable1to6-1))
                    .fill(time);
        }
        
    }
    
   
    private String[] parNames = {"par0", "par1"};
    private double[] errs = {0.1,0.1};
    
    
    public void loadFitPars() {
        
            double[] pars = new double[2];
            pars[0] = 2.0;
            pars[1] = 7.0;
            for (int i = 0; i < nsec; i++)
        {
            for (int j = 0; j < nsl; j++)
            {
                for (int k = 0; k < nSlots7; k++)
                {
                    for (int l = 0; l < nCables6; l++)
                    {
                        TDCFitPars.put(new Coordinate(i,j,k,l), new MnUserParameters());
                        for(int p = 0; p < parNames.length; p++) {
                            TDCFitPars.get(new Coordinate(i,j,k,l)).add(parNames[p], pars[p], errs[p]);
                        }
                    }
                }
            }
        }   
       
    }
    
    public void Plot(int i , int j, int k) {
        
        for (int l = 0; l < nCables6; l++){
        if(this.TDCHis.get(new Coordinate(i, j, k, l)).getEntries()>0) {
            this.getAnalysisCanvas().getCanvas(analTabs).cd(l);
            this.getAnalysisCanvas().getCanvas(analTabs)
                    .draw(this.TDCHis.get(new Coordinate(i, j, k, l)));
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

}

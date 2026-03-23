import javax.swing.*;
import javax.swing.border.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.geom.*;
import java.awt.image.BufferedImage;
import java.util.*;
import java.util.List;
import java.util.Timer;
import java.util.TimerTask;

public class GasSimulator extends JFrame {
    // Van der Waals constants for common gases
    static final double[][] GAS_CONSTANTS = {
        // {a (L²·atm/mol²), b (L/mol)}
        {3.640, 0.04267},  // CO2
        {1.360, 0.03183},  // N2
        {1.363, 0.03219},  // O2
        {2.253, 0.04281},  // Cl2
        {0.244, 0.02661},  // He
        {4.225, 0.03707},  // NH3
    };
    static final String[] GAS_NAMES = {"CO₂", "N₂", "O₂", "Cl₂", "He", "NH₃"};

    // Color scheme
    static final Color BG_DARK       = new Color(10, 12, 20);
    static final Color BG_PANEL      = new Color(18, 22, 36);
    static final Color BG_CARD       = new Color(24, 30, 50);
    static final Color ACCENT_CYAN   = new Color(0, 220, 255);
    static final Color ACCENT_ORANGE = new Color(255, 140, 0);
    static final Color ACCENT_PURPLE = new Color(160, 80, 255);
    static final Color ACCENT_GREEN  = new Color(60, 230, 120);
    static final Color TEXT_PRIMARY  = new Color(220, 230, 255);
    static final Color TEXT_DIM      = new Color(120, 140, 180);
    static final Color GRID_COLOR    = new Color(40, 50, 80);

    // UI Components
    private GraphPanel graphPanel;
    private SimulationPanel simPanel;
    private JComboBox<String> gasSelector;
    private JSlider tempSlider, molSlider, volSlider;
    private JLabel tempLabel, molLabel, volLabel;
    private JLabel zFactorLabel, pressureIdealLabel, pressureVdwLabel;
    private JPanel legendPanel;

    // Data
    private double temperature = 300;  // K
    private double moles = 1.0;        // mol
    private double volume = 1.0;       // L
    private int selectedGas = 0;

    public GasSimulator() {
        setTitle("Gas Behavior Simulator — Ideal vs. Van der Waals");
        setDefaultCloseOperation(EXIT_ON_CLOSE);
        setSize(1400, 860);
        setMinimumSize(new Dimension(1100, 700));
        setLocationRelativeTo(null);
        getContentPane().setBackground(BG_DARK);

        buildUI();
        setVisible(true);
        updateAll();
    }

    private void buildUI() {
        setLayout(new BorderLayout(8, 8));

        // Header 
        JPanel header = new JPanel(new BorderLayout());
        header.setBackground(BG_DARK);
        header.setBorder(new EmptyBorder(14, 20, 6, 20));

        JLabel title = new JLabel("IDEAL vs. NON-IDEAL GAS BEHAVIOR");
        title.setFont(new Font("Monospaced", Font.BOLD, 18));
        title.setForeground(ACCENT_CYAN);

        JLabel subtitle = new JLabel("Van der Waals Equation  ·  Compressibility Factor  ·  Molecular Simulation");
        subtitle.setFont(new Font("SansSerif", Font.PLAIN, 12));
        subtitle.setForeground(TEXT_DIM);

        JPanel titleBox = new JPanel(new GridLayout(2, 1, 0, 2));
        titleBox.setBackground(BG_DARK);
        titleBox.add(title);
        titleBox.add(subtitle);
        header.add(titleBox, BorderLayout.WEST);

        // Gas selector
        gasSelector = new JComboBox<>(GAS_NAMES);
        gasSelector.setBackground(BG_CARD);
        gasSelector.setForeground(TEXT_PRIMARY);
        gasSelector.setFont(new Font("Monospaced", Font.BOLD, 13));
        gasSelector.addActionListener(e -> { selectedGas = gasSelector.getSelectedIndex(); updateAll(); });
        JLabel gasLbl = new JLabel("Gas: ");
        gasLbl.setForeground(TEXT_DIM);
        gasLbl.setFont(new Font("SansSerif", Font.PLAIN, 12));
        JPanel gasBox = new JPanel(new FlowLayout(FlowLayout.RIGHT, 6, 0));
        gasBox.setBackground(BG_DARK);
        gasBox.add(gasLbl);
        gasBox.add(gasSelector);
        header.add(gasBox, BorderLayout.EAST);
        add(header, BorderLayout.NORTH);

        // Center => graph + simulation
        JSplitPane centerSplit = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT);
        centerSplit.setDividerSize(4);
        centerSplit.setBackground(BG_DARK);
        centerSplit.setDividerLocation(820);
        centerSplit.setBorder(null);

        graphPanel = new GraphPanel();
        centerSplit.setLeftComponent(wrap(graphPanel));

        simPanel = new SimulationPanel();
        centerSplit.setRightComponent(wrap(simPanel));
        add(centerSplit, BorderLayout.CENTER);

        // Bottom => controls + readouts
        JPanel bottom = new JPanel(new BorderLayout(10, 0));
        bottom.setBackground(BG_DARK);
        bottom.setBorder(new EmptyBorder(4, 12, 12, 12));

        bottom.add(buildControls(), BorderLayout.CENTER);
        bottom.add(buildReadouts(), BorderLayout.EAST);
        add(bottom, BorderLayout.SOUTH);
    }

    private JPanel wrap(JComponent c) {
        JPanel p = new JPanel(new BorderLayout());
        p.setBackground(BG_DARK);
        p.setBorder(BorderFactory.createLineBorder(new Color(40, 55, 90), 1));
        p.add(c);
        return p;
    }

    private JPanel buildControls() {
        JPanel p = new JPanel(new GridLayout(3, 1, 4, 4));
        p.setBackground(BG_DARK);
        p.setBorder(new EmptyBorder(4, 4, 4, 20));

        tempLabel = new JLabel();
        molLabel  = new JLabel();
        volLabel  = new JLabel();

        tempSlider = makeSlider(50, 1000, (int) temperature);
        molSlider  = makeSlider(1, 50, (int)(moles * 10));
        volSlider  = makeSlider(1, 200, (int)(volume * 10));

        tempSlider.addChangeListener(e -> { temperature = tempSlider.getValue(); updateAll(); });
        molSlider .addChangeListener(e -> { moles = molSlider.getValue() / 10.0;  updateAll(); });
        volSlider .addChangeListener(e -> { volume = volSlider.getValue() / 10.0; updateAll(); });

        p.add(makeLabeledSlider("Temperature (T)", tempLabel, tempSlider, "K"));
        p.add(makeLabeledSlider("Moles (n)", molLabel, molSlider, "mol"));
        p.add(makeLabeledSlider("Volume (V)", volLabel, volSlider, "L"));
        return p;
    }

    private JPanel makeLabeledSlider(String name, JLabel valLbl, JSlider slider, String unit) {
        JPanel row = new JPanel(new BorderLayout(8, 0));
        row.setBackground(BG_DARK);
        JLabel lbl = new JLabel(name);
        lbl.setFont(new Font("SansSerif", Font.PLAIN, 11));
        lbl.setForeground(TEXT_DIM);
        lbl.setPreferredSize(new Dimension(130, 20));
        valLbl.setFont(new Font("Monospaced", Font.BOLD, 12));
        valLbl.setForeground(ACCENT_CYAN);
        valLbl.setPreferredSize(new Dimension(80, 20));
        row.add(lbl, BorderLayout.WEST);
        row.add(slider, BorderLayout.CENTER);
        row.add(valLbl, BorderLayout.EAST);
        return row;
    }

    private JSlider makeSlider(int min, int max, int val) {
        JSlider s = new JSlider(min, max, val);
        s.setBackground(BG_DARK);
        s.setForeground(TEXT_DIM);
        return s;
    }

    private JPanel buildReadouts() {
        JPanel p = new JPanel(new GridLayout(3, 1, 4, 4));
        p.setBackground(BG_CARD);
        p.setBorder(new CompoundBorder(
            BorderFactory.createLineBorder(new Color(50, 65, 110), 1),
            new EmptyBorder(8, 14, 8, 14)));
        p.setPreferredSize(new Dimension(260, 90));

        pressureIdealLabel = makeReadout("P (Ideal):", ACCENT_CYAN);
        pressureVdwLabel   = makeReadout("P (VdW):", ACCENT_ORANGE);
        zFactorLabel       = makeReadout("Z (Compressibility):", ACCENT_GREEN);

        p.add(pressureIdealLabel);
        p.add(pressureVdwLabel);
        p.add(zFactorLabel);
        return p;
    }

    private JLabel makeReadout(String prefix, Color col) {
        JLabel l = new JLabel(prefix + " —");
        l.setFont(new Font("Monospaced", Font.BOLD, 12));
        l.setForeground(col);
        return l;
    }

    // Physics calculations
    static final double R = 0.082057; // L*atm / (mol*K)

    double idealPressure() {
        return (moles * R * temperature) / volume;
    }

    double vdwPressure() {
        double a = GAS_CONSTANTS[selectedGas][0];
        double b = GAS_CONSTANTS[selectedGas][1];
        double nb = moles * b;
        double denom = volume - nb;
        if (denom <= 0) return Double.NaN;
        return (moles * R * temperature / denom) - (a * moles * moles / (volume * volume));
    }

    double compressibilityZ() {
        double pi = idealPressure();
        double pv = vdwPressure();
        if (Double.isNaN(pv) || pi == 0) return Double.NaN;
        return pv / pi;
    }

    void updateAll() {
        double pi = idealPressure();
        double pv = vdwPressure();
        double z  = compressibilityZ();

        tempLabel.setText(String.format("%.0f K", temperature));
        molLabel .setText(String.format("%.1f mol", moles));
        volLabel .setText(String.format("%.1f L", volume));

        pressureIdealLabel.setText(String.format("P (Ideal):   %.4f atm", pi));
        pressureVdwLabel  .setText(String.format("P (VdW):     %.4f atm", Double.isNaN(pv) ? 0 : pv));
        zFactorLabel      .setText(String.format("Z factor:    %.4f", Double.isNaN(z) ? 0 : z));

        graphPanel.recompute(selectedGas, temperature, moles);
        simPanel.updateParams(temperature, moles, volume, z);
    }

    //  Graph Panel
    class GraphPanel extends JPanel {
        double[] volumes, pIdeal, pVdw, zFactors;
        int hoveredPoint = -1;
        static final int POINTS = 200;

        GraphPanel() {
            setBackground(BG_PANEL);
            setPreferredSize(new Dimension(820, 580));
            addMouseMotionListener(new MouseMotionAdapter() {
                public void mouseMoved(MouseEvent e) { findHover(e.getX(), e.getY()); }
            });
        }

        void recompute(int gas, double T, double n) {
            double a = GAS_CONSTANTS[gas][0];
            double b = GAS_CONSTANTS[gas][1];
            
            volumes  = new double[POINTS];
            pIdeal   = new double[POINTS];
            pVdw     = new double[POINTS];
            zFactors = new double[POINTS];
            
            double vMin = n * b * 1.05 + 0.01;
            double vMax = 30.0;
            
            for (int i = 0; i < POINTS; i++) {
                double v = vMin + (vMax - vMin) * i / (POINTS - 1);
                
                volumes[i] = v;
                pIdeal[i]  = n * R * T / v;
                
                double denom = v - n * b;
                
                pVdw[i] = denom > 0
                    ? (n * R * T / denom) - (a * n * n / (v * v))
                    : Double.NaN;
                    
                zFactors[i] = (Double.isNaN(pVdw[i]) || pIdeal[i] == 0)
                    ? Double.NaN
                    : pVdw[i] / pIdeal[i];
                    
            }
            repaint();
        }

        void findHover(int mx, int my) {
            if (volumes == null) return;
            Rectangle r = getPlotRect();
            if (!r.contains(mx, my)) { hoveredPoint = -1; repaint(); return; }
            double[] xRange = {volumes[0], volumes[POINTS-1]};
            double closest = Double.MAX_VALUE;
            int best = -1;
            for (int i = 0; i < POINTS; i++) {
                double px = mapX(volumes[i], xRange, r);
                double d  = Math.abs(px - mx);
                if (d < closest) { closest = d; best = i; }
            }
            hoveredPoint = (closest < 12) ? best : -1;
            repaint();
        }

        Rectangle getPlotRect() {
            int lm = 70, rm = 30, tm = 50, bm = 60;
            return new Rectangle(lm, tm, getWidth() - lm - rm, getHeight() / 2 - tm - bm / 2);
        }
        Rectangle getZRect() {
            int lm = 70, rm = 30, bm = 50;
            int h = getHeight() / 2;
            return new Rectangle(lm, h + 20, getWidth() - lm - rm, h - 50);
        }

        double mapX(double v, double[] xr, Rectangle r) {
            return r.x + (v - xr[0]) / (xr[1] - xr[0]) * r.width;
        }
        double mapY(double p, double[] yr, Rectangle r) {
            return r.y + r.height - (p - yr[0]) / (yr[1] - yr[0]) * r.height;
        }

        @Override
        protected void paintComponent(Graphics g0) {
            super.paintComponent(g0);
            Graphics2D g = (Graphics2D) g0;
            g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

            drawBackground(g);
            if (volumes == null) return;

            Rectangle pr = getPlotRect();
            Rectangle zr = getZRect();

            // Pressure vs Volume chart
            drawGrid(g, pr, "Pressure (atm)", "Volume (L)");
            double[] xr = {volumes[0], volumes[POINTS-1]};
            double pMax = 0;
            for (double p : pIdeal) if (!Double.isNaN(p) && p < 100) pMax = Math.max(pMax, p);
            double[] yr = {0, pMax * 1.1};

            drawCurve(g, volumes, pIdeal, xr, yr, pr, ACCENT_CYAN, new float[]{6,0}, "Ideal Gas Law");
            drawCurve(g, volumes, pVdw,   xr, yr, pr, ACCENT_ORANGE, new float[]{8,4}, "Van der Waals");
            drawAxisLabels(g, pr, yr, "Pressure (atm)");
            drawChartTitle(g, pr, "Pressure vs. Volume  (P-V Diagram)");

            // Z factor chart
            drawGrid(g, zr, "Compressibility Z = P_VdW / P_Ideal", "Volume (L)");
            double[] zr2 = {0.5, 1.5};
            // draw Z=1 reference line
            g.setColor(new Color(255, 255, 255, 50));
            g.setStroke(new BasicStroke(1.2f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND,
                0, new float[]{4,4}, 0));
            int y1 = (int) mapY(1.0, zr2, zr);
            g.drawLine(zr.x, y1, zr.x + zr.width, y1);
            g.setColor(new Color(255,255,255,80));
            g.setFont(new Font("Monospaced", Font.PLAIN, 10));
            g.drawString("Z = 1  (ideal)", zr.x + 4, y1 - 4);

            drawCurve(g, volumes, zFactors, xr, zr2, zr, ACCENT_GREEN, new float[]{6,0}, "Z Factor");
            drawAxisLabels(g, zr, zr2, "Compressibility Z");
            drawChartTitle(g, zr, "Compressibility Factor Z  (Deviation from Ideal)");

            // Legend
            drawLegend(g, zr);

            // Hover tooltip
            if (hoveredPoint >= 0) drawTooltip(g, hoveredPoint, xr, yr, pr);
        }

        void drawBackground(Graphics2D g) {
            g.setColor(BG_PANEL);
            g.fillRect(0, 0, getWidth(), getHeight());
            // subtle gradient at top
            GradientPaint gp = new GradientPaint(0, 0, new Color(20, 35, 70, 80),
                0, getHeight() / 2, new Color(0, 0, 0, 0));
            g.setPaint(gp);
            g.fillRect(0, 0, getWidth(), getHeight());
        }

        void drawGrid(Graphics2D g, Rectangle r, String title, String xLabel) {
            g.setColor(BG_CARD);
            g.fillRect(r.x, r.y, r.width, r.height);
            g.setColor(GRID_COLOR);
            g.setStroke(new BasicStroke(0.5f));
            for (int i = 1; i <= 8; i++) {
                int x = r.x + r.width * i / 8;
                int y = r.y + r.height * i / 8;
                g.drawLine(x, r.y, x, r.y + r.height);
                g.drawLine(r.x, y, r.x + r.width, y);
            }
            g.setColor(new Color(60, 80, 120));
            g.setStroke(new BasicStroke(1));
            g.drawRect(r.x, r.y, r.width, r.height);
        }

        void drawCurve(Graphics2D g, double[] xs, double[] ys, double[] xr, double[] yr,
                        Rectangle r, Color col, float[] dash, String label) {
            g.setColor(col);
            g.setStroke(new BasicStroke(2.2f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND,
                0, dash, 0));
            GeneralPath path = new GeneralPath();
            boolean started = false;
            for (int i = 0; i < xs.length; i++) {
                if (Double.isNaN(ys[i])) { started = false; continue; }
                double px = mapX(xs[i], xr, r);
                double py = mapY(ys[i], yr, r);
                
                if (py < r.y - 10 || py > r.y + r.height + 10) { 
                    started = false; 
                    continue; 
                }
                
                if (!started) { 
                    path.moveTo(px, py); 
                    started = true; 
                }
                else path.lineTo(px, py);
            }
            g.draw(path);
        }

        void drawAxisLabels(Graphics2D g, Rectangle r, double[] yr, String yLabel) {
            g.setFont(new Font("Monospaced", Font.PLAIN, 10));
            g.setColor(TEXT_DIM);
            for (int i = 0; i <= 4; i++) {
                double val = yr[0] + (yr[1] - yr[0]) * i / 4;
                int y = (int)(r.y + r.height - r.height * i / 4);
                g.drawString(String.format("%.2f", val), r.x - 52, y + 4);
                g.setColor(GRID_COLOR);
            }
            // X axis labels
            g.setColor(TEXT_DIM);
            double vMin = volumes[0], vMax = volumes[POINTS - 1];
            for (int i = 0; i <= 5; i++) {
                double v = vMin + (vMax - vMin) * i / 5;
                int x = (int) mapX(v, new double[]{vMin, vMax}, r);
                g.drawString(String.format("%.1f", v), x - 10, r.y + r.height + 16);
            }
            // Rotated Y label
            Graphics2D g2 = (Graphics2D) g.create();
            g2.setColor(TEXT_DIM);
            g2.setFont(new Font("SansSerif", Font.PLAIN, 11));
            g2.translate(r.x - 62, r.y + r.height / 2);
            g2.rotate(-Math.PI / 2);
            g2.drawString(yLabel, -yLabel.length() * 3, 0);
            g2.dispose();

            g.setColor(TEXT_DIM);
            g.setFont(new Font("SansSerif", Font.PLAIN, 11));
            g.drawString("Volume (L)", r.x + r.width / 2 - 28, r.y + r.height + 32);
        }

        void drawChartTitle(Graphics2D g, Rectangle r, String title) {
            g.setFont(new Font("SansSerif", Font.BOLD, 12));
            g.setColor(TEXT_PRIMARY);
            g.drawString(title, r.x + 10, r.y - 10);
        }

        void drawLegend(Graphics2D g, Rectangle r) {
            int lx = r.x + r.width - 250, ly = r.y + r.height + 10;
            String[][] items = {
                {"Ideal Gas  (PV = nRT)", String.valueOf(ACCENT_CYAN.getRGB())},
                {"Van der Waals  (non-ideal)", String.valueOf(ACCENT_ORANGE.getRGB())},
                {"Compressibility Z", String.valueOf(ACCENT_GREEN.getRGB())},
            };
            g.setFont(new Font("SansSerif", Font.PLAIN, 11));
            for (int i = 0; i < items.length; i++) {
                Color c = new Color(Integer.parseInt(items[i][1]));
                g.setColor(c);
                g.setStroke(new BasicStroke(2.5f));
                g.drawLine(lx + i * 230 - 440, ly + 2, lx + i * 230 - 415, ly + 2);
                g.setColor(TEXT_DIM);
                g.drawString(items[i][0], lx + i * 230 - 410, ly + 6);
            }
        }

        void drawTooltip(Graphics2D g, int idx, double[] xr, double[] yr, Rectangle r) {
            double v  = volumes[idx];
            double pi = pIdeal[idx];
            double pv = Double.isNaN(pVdw[idx]) ? 0 : pVdw[idx];
            double z  = Double.isNaN(zFactors[idx]) ? 0 : zFactors[idx];
            int tx = (int) mapX(v, xr, r) + 10;
            int ty = r.y + 20;

            String[] lines = {
                String.format("V = %.3f L", v),
                String.format("P_ideal = %.4f atm", pi),
                String.format("P_VdW   = %.4f atm", pv),
                String.format("Z       = %.4f", z),
            };
            int tw = 175, th = lines.length * 16 + 14;
            if (tx + tw > r.x + r.width) tx = (int) mapX(v, xr, r) - tw - 10;

            g.setColor(new Color(10, 15, 30, 210));
            g.fillRoundRect(tx, ty, tw, th, 8, 8);
            g.setColor(ACCENT_CYAN);
            g.setStroke(new BasicStroke(1));
            g.drawRoundRect(tx, ty, tw, th, 8, 8);
            g.setFont(new Font("Monospaced", Font.PLAIN, 11));
            for (int i = 0; i < lines.length; i++) {
                g.setColor(i == 0 ? ACCENT_CYAN : TEXT_PRIMARY);
                g.drawString(lines[i], tx + 8, ty + 14 + i * 16);
            }
        }
    }

    //  Simulation Panel -> animated particles
    class SimulationPanel extends JPanel {
        private List<Particle> particles = new ArrayList<>();
        private Timer animTimer;
        private double temp = 300, n = 1.0, vol = 1.0, z = 1.0;
        private int containerH = 300;
        private Random rng = new Random();
        private int tick = 0;

        SimulationPanel() {
            setBackground(BG_PANEL);
            setPreferredSize(new Dimension(380, 580));
            initParticles(20);
            animTimer = new Timer(true);
            animTimer.scheduleAtFixedRate(new TimerTask() {
                public void run() {
                    tick++;
                    for (Particle p : particles) p.update();
                    SwingUtilities.invokeLater(() -> repaint());
                }
            }, 0, 33);
        }

        void updateParams(double t, double nm, double v, double zz) {
            temp = t; n = nm; vol = v; z = Double.isNaN(zz) ? 1.0 : zz;
            int count = Math.min(5 + (int)(nm * 4), 60);
            initParticles(count);
        }

        void initParticles(int count) {
            particles.clear();
            for (int i = 0; i < count; i++) particles.add(new Particle(rng));
        }

        @Override
        protected void paintComponent(Graphics g0) {
            super.paintComponent(g0);
            Graphics2D g = (Graphics2D) g0;
            g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

            int W = getWidth(), H = getHeight();
            g.setColor(BG_PANEL);
            g.fillRect(0, 0, W, H);

            // Title
            g.setFont(new Font("SansSerif", Font.BOLD, 13));
            g.setColor(TEXT_PRIMARY);
            g.drawString("Molecular Simulation", 16, 24);
            g.setFont(new Font("SansSerif", Font.PLAIN, 11));
            g.setColor(TEXT_DIM);
            g.drawString("Particle speeds ∝ temperature  ·  Size = molecular volume", 16, 40);

            // Container
            int cx = 20, cy = 55, cw = W - 40;
            containerH = Math.max(150, Math.min(H - 200, (int)(50 + vol * 18)));
            g.setColor(new Color(30, 45, 80));
            g.fillRect(cx, cy, cw, containerH);
            g.setColor(new Color(60, 90, 160));
            g.setStroke(new BasicStroke(2));
            g.drawRect(cx, cy, cw, containerH);

            // Pressure glow on walls based on Z
            float alpha = (float) Math.min(0.6, Math.abs(z - 1.0) * 1.5);
            Color wallColor = z > 1.0 ? new Color(255, 140, 0, (int)(alpha * 255))
                                       : new Color(0, 220, 255, (int)(alpha * 255));
            g.setColor(wallColor);
            g.setStroke(new BasicStroke(3));
            g.drawRect(cx, cy, cw, containerH);

            // Particles
            for (Particle p : particles) {
                p.clamp(cx + 4, cy + 4, cx + cw - 4, cy + containerH - 4);
                p.draw(g, temp, z);
            }

            // Legend box
            drawColorLegend(g, cx, cy + containerH + 14, cw);

            // Stats
            drawStats(g, cx, cy + containerH + 100);
        }

        void drawColorLegend(Graphics2D g, int x, int y, int w) {
            g.setColor(BG_CARD);
            g.fillRoundRect(x, y, w, 70, 8, 8);
            g.setColor(new Color(50, 65, 110));
            g.setStroke(new BasicStroke(1));
            g.drawRoundRect(x, y, w, 70, 8, 8);

            g.setFont(new Font("SansSerif", Font.BOLD, 11));
            g.setColor(TEXT_PRIMARY);
            g.drawString("Particle Color Key:", x + 10, y + 16);

            String[][] keys = {
                {"Core color", "Temperature (blue=cold → red=hot)"},
                {"Ring color", "Kinetic energy (green=low → yellow=high)"},
                {"Glow", "Intermolecular interaction strength"},
            };
            g.setFont(new Font("SansSerif", Font.PLAIN, 10));
            for (int i = 0; i < keys.length; i++) {
                g.setColor(TEXT_DIM);
                g.drawString("● " + keys[i][0] + ":  " + keys[i][1], x + 10, y + 30 + i * 14);
            }
        }

        void drawStats(Graphics2D g, int x, int y) {
            g.setFont(new Font("Monospaced", Font.PLAIN, 11));
            String[] stats = {
                String.format("Particles: %d", particles.size()),
                String.format("T = %.0f K   n = %.1f mol   V = %.1f L", temp, n, vol),
                String.format("Z = %.4f  (%s)", z,
                    z > 1.01 ? "repulsive dominates" : z < 0.99 ? "attractive dominates" : "near-ideal"),
            };
            Color[] colors = {ACCENT_CYAN, TEXT_PRIMARY, z > 1.01 ? ACCENT_ORANGE : z < 0.99 ? ACCENT_CYAN : ACCENT_GREEN};
            for (int i = 0; i < stats.length; i++) {
                g.setColor(colors[i]);
                g.drawString(stats[i], x, y + i * 16);
            }
        }

        class Particle {
            double x, y, vx, vy;
            double baseSpeed;
            double phase;
            int size;

            Particle(Random rng) {
                x = 30 + rng.nextDouble() * 300;
                y = 60 + rng.nextDouble() * 200;
                phase = rng.nextDouble() * Math.PI * 2;
                size = 5 + rng.nextInt(5);
                double angle = rng.nextDouble() * Math.PI * 2;
                baseSpeed = 0.8 + rng.nextDouble() * 1.5;
                vx = Math.cos(angle) * baseSpeed;
                vy = Math.sin(angle) * baseSpeed;
            }

            void update() {
                double speedScale = Math.sqrt(temp / 300.0);
                x += vx * speedScale;
                y += vy * speedScale;
            }

            void clamp(int x1, int y1, int x2, int y2) {
                if (x < x1) { x = x1; vx = Math.abs(vx); }
                if (x > x2) { x = x2; vx = -Math.abs(vx); }
                if (y < y1) { y = y1; vy = Math.abs(vy); }
                if (y > y2) { y = y2; vy = -Math.abs(vy); }
            }

            void draw(Graphics2D g, double temp, double z) {
                // Core color: temperature (blue → red)
                float t = (float) Math.min(1, Math.max(0, (temp - 50) / 950.0));
                Color core = new Color(
                    (int)(t * 255),
                    (int)(Math.sin(t * Math.PI) * 180),
                    (int)((1 - t) * 255)
                );

                // Ring color: kinetic energy proxy
                double speed = Math.sqrt(vx * vx + vy * vy) * Math.sqrt(temp / 300.0);
                float e = (float) Math.min(1, speed / 4.0);
                Color ring = new Color(
                    (int)(e * 255),
                    (int)((1 - e * 0.5) * 230),
                    0
                );

                // Glow: intermolecular forces (from Z deviation)
                double zDev = Math.abs(z - 1.0);
                int glowAlpha = (int) Math.min(180, zDev * 400);
                Color glowColor = z > 1.0
                    ? new Color(255, 100, 0, glowAlpha)
                    : new Color(0, 180, 255, glowAlpha);

                int xi = (int) x, yi = (int) y;
                int r = size;

                // Glow halo
                if (glowAlpha > 10) {
                    for (int gw = r + 8; gw > r; gw--) {
                        int ga = glowAlpha * (gw - r) / 8;
                        g.setColor(new Color(glowColor.getRed(), glowColor.getGreen(),
                            glowColor.getBlue(), Math.max(0, 50 - ga)));
                        g.fillOval(xi - gw, yi - gw, gw * 2, gw * 2);
                    }
                }

                // Core body
                RadialGradientPaint rgp = new RadialGradientPaint(
                    new Point2D.Float(xi - r / 3f, yi - r / 3f),
                    r,
                    new float[]{0f, 1f},
                    new Color[]{core.brighter(), core.darker()}
                );
                g.setPaint(rgp);
                g.fillOval(xi - r, yi - r, r * 2, r * 2);

                // Ring
                g.setColor(ring);
                g.setStroke(new BasicStroke(1.5f));
                g.drawOval(xi - r, yi - r, r * 2, r * 2);

                // Highlight
                g.setColor(new Color(255, 255, 255, 80));
                g.fillOval(xi - r / 3, yi - r / 2, r / 2, r / 3);
            }
        }
    }

    public static void main(String[] args) {
        System.setProperty("awt.useSystemAAFontSettings", "on");
        System.setProperty("swing.aatext", "true");
        try {
            UIManager.setLookAndFeel(
                UIManager.getSystemLookAndFeelClassName()
            );
        }
        catch (Exception ignored) {
            
        }
        SwingUtilities.invokeLater(GasSimulator::new);
    }
}

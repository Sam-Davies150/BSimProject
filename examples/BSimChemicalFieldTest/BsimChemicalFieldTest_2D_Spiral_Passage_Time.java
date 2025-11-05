package BSimChemicalFieldTest;

import bsim.BSim;
import bsim.BSimChemicalField;
import bsim.BSimTicker;
import bsim.BSimUtils;
import bsim.draw.BSimP3DDrawer;
import bsim.export.BSimLogger;
import bsim.particle.BSimBacterium;
import processing.core.PGraphics3D;

import javax.vecmath.Vector3d;
import java.awt.*;
import java.util.ArrayList;
import java.util.Vector;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.List;
import flanagan.integration.Integration;
import flanagan.integration.IntegralFunction;


public class BsimChemicalFieldTest_2D_Spiral_Passage_Time {

    // Helper class to store spiral position information
    static class SpiralPoint {
        double angle;
        Vector3d position;

        SpiralPoint(double angle, Vector3d position) {
            this.angle = angle;
            this.position = position;
        }
    }

    // Class to track passage time v2 - pls work
    static class PassageTimeTracker {
        private Double startBacteriumFireTime = null;
        private Double endBacteriumFireTime = null;
        private Double passageTime = null;

        synchronized void setStartFire(double time) {
            if (this.startBacteriumFireTime == null) {
                this.startBacteriumFireTime = time;
                checkPassageTimeComplete();
            }
        }

        synchronized void setEndFire(double time) {
            if (this.endBacteriumFireTime == null) {
                this.endBacteriumFireTime = time;
                checkPassageTimeComplete();
            }
        }

        private void checkPassageTimeComplete() {
            if (this.startBacteriumFireTime != null && this.endBacteriumFireTime != null) {
                this.passageTime = this.endBacteriumFireTime - this.startBacteriumFireTime;
                System.out.println("PASSAGE TIME CALCULATED: " + this.passageTime + " seconds");
            }
        }

        synchronized Double getPassageTime() {
            return this.passageTime;
        }
    }

    // Adding a class below to implement the arc length calculation for the integral
    static class SpiralArcLengthFunction implements IntegralFunction {
        private double a;
        public SpiralArcLengthFunction(double a) {
            this.a = a;
        }

        @Override
        public double function(double theta) {
            // f(θ) = sqrt((aθ)^2 + a^2)
            return Math.sqrt(Math.pow(a * theta, 2) + Math.pow(a, 2));
        }
    }

    public static void run(double decayRate,double diffusivity,double threshold,double quantityadd,int subdivisionLevels) {
        final int bounds = 110;
        BSim sim = new BSim();
        sim.setDt(0.002);
        sim.setTimeFormat("0.000");
        sim.setSolid(true, true, true);
        sim.setBound(bounds, bounds, bounds);
        if(subdivisionLevels > 9){sim.setSimulationTime(10);}
        else if(subdivisionLevels > 5){sim.setSimulationTime(50);}
        else {sim.setSimulationTime(100);}

        final int boxcount = 32;

//        final double decayRate = 0.00044;
//        final double diffusivity = 0.5;
//        final double threshold = 3e6;
//        final double quantityadd = 1e12;
        final double firstwave = 180;
        final double secondwave = 1440;

        final BSimChemicalField field = new BSimChemicalField(sim, new int[]{boxcount, boxcount, boxcount}, diffusivity, decayRate);

        final Vector<BSimBacterium> bacteria = new Vector<>();

        // Archimedes Spiral parameters for 2D placement
        final double centerX = bounds / 2.0;
        final double centerY = bounds / 2.0;
        final double constantZ = bounds / 2.0;
        final double maxRadius = bounds * 0.4;
        final double spiralTurns = 1;

        // Archimedes spiral formula: r = a + b*theta
        final double a = 2.0;
        final double b = maxRadius / (2 * Math.PI * spiralTurns);
        final double spiralThickness = 5.0;

        // Parameters for calculating the volume of the spiral - gonna do it as a coiled cylinder
        final double radius = spiralThickness/2;

        Random rand = new Random();

        // NEW: Subdivision parameters
//        final int subdivisionLevels = 5; // Number of recursive subdivisions (adjust for more/fewer bacteria)
        // subdivisionLevels = 1 -> 2 bacteria (start, end)
        // subdivisionLevels = 2 -> 3 bacteria (start, middle, end)
        // subdivisionLevels = 3 -> 5 bacteria
        // subdivisionLevels = 4 -> 9 bacteria
        // subdivisionLevels = 5 -> 17 bacteria
        // Each level adds 2^(level-1) new bacteria

        // Define start and end angles for the spiral
        final double startAngle = 0.0;
        final double endAngle = 2 * Math.PI * spiralTurns;

        // Create list to store spiral points
        List<SpiralPoint> spiralPoints = new ArrayList<>();

        // Create passage time tracker
        final PassageTimeTracker passageTimeTracker = new PassageTimeTracker();

        // Recursive function to add bacteria at midpoints
        addBacteriaRecursive(spiralPoints, startAngle, endAngle, subdivisionLevels,
                subdivisionLevels, a, b, centerX, centerY, constantZ, spiralThickness,
                bounds, rand, sim, bacteria, field, quantityadd,
                firstwave, secondwave, threshold, passageTimeTracker,
                startAngle, endAngle);

        System.out.println("Total bacteria created: " + bacteria.size());

        // Below is the actual integration for the arc length of the spiral for the volume calc

        double theta1 = 0.0;
        double theta2 = 2 * Math.PI;

        SpiralArcLengthFunction func = new SpiralArcLengthFunction(a);

// Create an Integration object (using Gauss-Legendre with 100 points)
        Integration integral = new Integration(func, theta1, theta2);
        double L = integral.gaussQuad(100);
// Calculate spiral volume
        double spiral_volume = Math.PI * Math.pow(radius,2) * L;
        System.out.println("Spiral Volume: " + spiral_volume);
// Calculate number density (pls work)
        double number_density = bacteria.size()/spiral_volume;
        System.out.println("Number Density: " + number_density);


        final int numThreads = Runtime.getRuntime().availableProcessors();
        final ExecutorService executorService = Executors.newFixedThreadPool(numThreads);

        sim.setTicker(new BSimTicker() {
            public void tick() {
                bacteria.parallelStream().forEach(b -> b.action());
                field.update();
            }
        });

        sim.setDrawer(new BSimP3DDrawer(sim, 800, 800) {
            @Override
            public void scene(PGraphics3D p3d) {
                draw(field, new Color(0, 255, 255), (float) (10000 / quantityadd));
                for (BSimBacterium p : bacteria) draw(p, Color.RED);
            }
        });

        String resultsDir = BSimUtils.generateDirectoryPath("./results/");
        String resultsPath = resultsDir + "Passage_Time4.csv";
        BSimLogger logger = new BSimLogger(sim, resultsPath) {
            @Override
            public void before() {
                super.before();
                java.io.File f = new java.io.File(resultsPath);
                if (f.length() == 0) {
                    StringBuilder header = new StringBuilder("decayRate,diffusivity,threshold,quantityAdd,subdivisionLevel,Number,PassageTime");
                    write(header.toString());
                }
            }


            @Override
            public void during() {
                if(sim.getTime() == sim.getSimulationTime()) {
                    StringBuilder line = new StringBuilder();
                    line.append(decayRate).append(",").
                            append(diffusivity).append(",").
                            append(threshold).append(",").
                            append(quantityadd).append(",").
                            append(subdivisionLevels).append(",").
                            append(bacteria.size()).append(",")
                            .append(passageTimeTracker.getPassageTime());

                    write(line.toString());
                }
            }
        };

        sim.addExporter(logger);

        Runtime.getRuntime().addShutdownHook(new Thread(() -> {
            executorService.shutdown();
            try {
                if (!executorService.awaitTermination(5, TimeUnit.SECONDS)) {
                    executorService.shutdownNow();
                }
            } catch (InterruptedException e) {
                executorService.shutdownNow();
            }
        }));

        //sim.preview();
        sim.export();
    }

    // Recursive function to create bacteria at midpoints along the spiral
    private static void addBacteriaRecursive(List<SpiralPoint> spiralPoints,
                                             double angleStart, double angleEnd, int depth,
                                             int maxDepth,
                                             double a, double b, double centerX, double centerY,
                                             double constantZ, double spiralThickness, int bounds,
                                             Random rand, BSim sim, Vector<BSimBacterium> bacteria,
                                             BSimChemicalField field, double quantityadd,
                                             double firstwave, double secondwave, double threshold,
                                             PassageTimeTracker tracker, double globalStart, double globalEnd) {
        if (depth == 0) {
            return;
        }

        // Create bacteria at start and end points
        if (depth == maxDepth) {
            createBacteriumAtAngle(angleStart, a, b, centerX, centerY, constantZ,
                    spiralThickness, bounds, rand, sim, bacteria,
                    field, quantityadd, firstwave, secondwave, threshold,
                    tracker, true, false);
            createBacteriumAtAngle(angleEnd, a, b, centerX, centerY, constantZ,
                    spiralThickness, bounds, rand, sim, bacteria,
                    field, quantityadd, firstwave, secondwave, threshold,
                    tracker, false, true);
        }
        if (depth > 1) {
            // Calculate midpoint angle
            double midAngle = (angleStart + angleEnd) / 2.0;

            // Create bacterium at midpoint
            createBacteriumAtAngle(midAngle, a, b, centerX, centerY, constantZ,
                    spiralThickness, bounds, rand, sim, bacteria,
                    field, quantityadd, firstwave, secondwave, threshold,
                    tracker, false, false);

            // Recursively subdivide the two halves
            addBacteriaRecursive(spiralPoints, angleStart, midAngle, depth - 1,
                    maxDepth, a, b, centerX, centerY, constantZ, spiralThickness,
                    bounds, rand, sim, bacteria, field, quantityadd,
                    firstwave, secondwave, threshold, tracker, globalStart, globalEnd);

            addBacteriaRecursive(spiralPoints, midAngle, angleEnd, depth - 1,
                    maxDepth, a, b, centerX, centerY, constantZ, spiralThickness,
                    bounds, rand, sim, bacteria, field, quantityadd,
                    firstwave, secondwave, threshold, tracker, globalStart, globalEnd);
        }
    }

    // Helper function to create a bacterium at a specific angle on the spiral
    private static void createBacteriumAtAngle(double angle, double a, double b,
                                               double centerX, double centerY, double constantZ,
                                               double spiralThickness, int bounds, Random rand,
                                               BSim sim, Vector<BSimBacterium> bacteria,
                                               BSimChemicalField field, double quantityadd,
                                               double firstwave, double secondwave, double threshold,
                                               PassageTimeTracker tracker, boolean isStart, boolean isEnd) {
        // Calculate radius at this angle
        double radius = a + b * angle;

        // Add random offset within spiral thickness
        double radialOffset = (rand.nextDouble() - 0.5) * spiralThickness;
        double effectiveRadius = radius + radialOffset;

        // Small angular variation for natural distribution
        double angularOffset = (rand.nextDouble() - 0.5) * 0.1;
        double effectiveAngle = angle;// + angularOffset;

        double x = centerX + effectiveRadius * Math.cos(effectiveAngle);
        double y = centerY + effectiveRadius * Math.sin(effectiveAngle);
        double z = constantZ + (rand.nextDouble() - 0.5) * spiralThickness;

        // Check bounds
        if (x < 0 || x >= bounds || y < 0 || y >= bounds) {
            return;
        }

        BSimBacterium p = new BSimBacterium(sim, new Vector3d(x, y, z)) {
            public boolean isFiring = false;
            public double fireStartTime = -1;
            public double fireDuration = 0.1;
            public boolean hasFiredFirstWave = false;

            public boolean isFiringSecondWave = false;
            public double fireStartTime2 = -1;
            public double fireDuration2 = 0.1;
            public boolean hasFiredSecondWave = false;

            @Override
            public void action() {
                super.action();
                double t = sim.getTime();

                double fireRate = quantityadd / (fireDuration / sim.getDt());
                double fireRate2 = quantityadd / (fireDuration2 / sim.getDt());

                // Initial spike in center
                if (t < 7 * sim.getDt()) {
                    Vector3d center = new Vector3d(bounds / 2.0, bounds / 2.0, bounds / 2.0);
                    double dx = position.x - center.x;
                    double dy = position.y - center.y;
                    double dz = position.z - center.z;
                    double dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
                    if (dist < 5.0) {
                        field.addQuantity(position, quantityadd);
                    }
                }

                // First wave firing logic
                if (!hasFiredFirstWave &&
                        t < firstwave &&
                        field.getConc(position) > threshold) {

                    fireStartTime = t;
                    isFiring = true;
                    hasFiredFirstWave = true;
                    setHasFired(1);

                    // Track passage time
                    if (isStart) {
                        tracker.setStartFire(t);
                    }
                    if (isEnd) {
                        tracker.setEndFire(t);
                    }
                }

                if (isFiring && (t - fireStartTime) <= fireDuration) {
                    field.addQuantity(position, quantityadd);
                } else {
                    isFiring = false;
                }

                // Second wave firing logic
                if (!hasFiredSecondWave &&
                        t > secondwave &&
                        field.getConc(position) > threshold) {

                    fireStartTime2 = t;
                    isFiringSecondWave = true;
                    hasFiredSecondWave = true;
                }

                if (isFiringSecondWave && (t - fireStartTime2) <= fireDuration2) {
                    field.addQuantity(position, quantityadd);
                } else {
                    isFiringSecondWave = false;
                }
            }
        };

        p.setGoal(field);
        if (!p.intersection(bacteria)) {
            bacteria.add(p);
        }
    }
    public static void main(String[] args){
        double[] decayRates = {0.00044};
        double[] diffusivities = {0.5};
        double[] thresholds = {3e6};
        double[] quantityAdds = {1e12};
        int MaxSubdivisions = 10;


        for (double decayRate : decayRates) {
            for (double diffusivity : diffusivities) {
                for (double threshold : thresholds) {
                    for (double quantityAdd : quantityAdds) {
                        for (int level = 1; level <= MaxSubdivisions; level++) {
                            System.out.println("Running simulation with:");
                            System.out.println("decayRate=" + decayRate +
                                    ", diffusivity=" + diffusivity +
                                    ", threshold=" + threshold +
                                    ", quantityAdd=" + quantityAdd +
                                    ", level=" + level);
                            run(decayRate, diffusivity, threshold, quantityAdd, level);
                        }
                    }
                }
            }
        }
    }
}
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
import java.util.Vector;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.IntStream;

public class BSimChemicalFieldTest_2D_Spiral {

    public static void main(String[] args) {
        final int bounds = 110;
        BSim sim = new BSim();
        sim.setDt(0.002);
        sim.setTimeFormat("0.000");
        sim.setSolid(true, true, true);
        sim.setBound(bounds, bounds, bounds);
        sim.setSimulationTime(120);

        final int boxcount = 64;

        final double decayRate = 0.07;
        final double diffusivity = 1380;
        final double threshold = 3e8;
        final double quantityadd = 3e10;
        final double firstwave = 180;
        final double secondwave = 1440;

        final BSimChemicalField field = new BSimChemicalField(sim, new int[]{boxcount, boxcount, boxcount}, diffusivity, decayRate);

        final Vector<BSimBacterium> bacteria = new Vector<>();
        int totbac = 100;

        // Archimedes Spiral parameters for 2D placement
        final double centerX = bounds / 2.0;
        final double centerY = bounds / 2.0;
        final double constantZ = bounds / 2.0;
        final double maxRadius = bounds * 0.4;
        final double spiralTurns = 1;

        // Archimedes spiral formula: r = a + b*theta
        final double a = 2.0;
        final double b = maxRadius / (2 * Math.PI * spiralTurns);

        // ADDED: Spiral thickness parameter (adjust this value to change thickness)
        final double spiralThickness = 5.0; // Thickness in microns (radial spread)

        Random rand = new Random();
        double angle = 0;
        double angleIncrement = (2 * Math.PI * spiralTurns) / totbac;

        while (bacteria.size() < totbac) {
            // Archimedes spiral coordinates: r = a + b*theta
            double radius = a + b * angle;

            if (radius > maxRadius) {
                angle = 0;
                angleIncrement *= 0.9;
                continue;
            }

            // ADDED: Random offset within spiral thickness
            // Random radial offset from spiral centerline
            double radialOffset = (rand.nextDouble() - 0.5) * spiralThickness;
            double effectiveRadius = radius + radialOffset;

            // Random angular perturbation (optional, for more natural distribution)
            double angularOffset = (rand.nextDouble() - 0.5) * 0.1; // Small angular variation
            double effectiveAngle = angle + angularOffset;

            double x = centerX + effectiveRadius * Math.cos(effectiveAngle);
            double y = centerY + effectiveRadius * Math.sin(effectiveAngle);
            double z = constantZ;

            // Check bounds to prevent out of bounds bacteria
            if (x < 0 || x >= bounds || y < 0 || y >= bounds) {
                angle += angleIncrement;
                continue;
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

            angle += angleIncrement;
        }

        final int numThreads = Runtime.getRuntime().availableProcessors();
        final ExecutorService executorService = Executors.newFixedThreadPool(numThreads);
// Paralellised ticker below
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
        BSimLogger logger = new BSimLogger(sim, resultsDir + "Archimedes_Spiral_Intensity.csv") {
            @Override
            public void before() {
                super.before();
                StringBuilder header = new StringBuilder("time");
                int maxShells = boxcount / 2;
                for (int r = 0; r < maxShells; r++) header.append(",r").append(r);
                write(header.toString());
            }

            @Override
            public void during() {
                double voxelSize_um = sim.getBound().x / boxcount;
                double centerX_um = sim.getBound().x / 2.0;
                double centerY_um = sim.getBound().y / 2.0;
                double centerZ_um = sim.getBound().z / 2.0;
                int maxShells = boxcount / 2;

                StringBuilder line = new StringBuilder(sim.getFormattedTime());
// Below is parallelised using claude
                double[] shellAverages = IntStream.range(0, maxShells)
                        .parallel()
                        .mapToDouble(r_vox -> {
                            double r_min = r_vox * voxelSize_um;
                            double r_max = (r_vox + 1) * voxelSize_um;
                            double sum = 0.0;
                            int count = 0;

                            for (int x = 0; x < boxcount; x++) {
                                for (int y = 0; y < boxcount; y++) {
                                    int z = 0;
                                    double x_um = (x + 0.5) * voxelSize_um;
                                    double y_um = (y + 0.5) * voxelSize_um;
                                    double z_um = (z + 0.5) * voxelSize_um;
                                    double dx = x_um - centerX_um;
                                    double dy = y_um - centerY_um;
                                    double r = Math.sqrt(dx * dx + dy * dy);

                                    if (r >= r_min && r < r_max) {
                                        double c = field.getConc(x, y, z);
                                        if (Double.isFinite(c) && c < 1e20) {
                                            sum += c;
                                            count++;
                                        }
                                    }
                                }
                            }

                            double avg = (count > 0) ? sum / count : 0.0;
                            return Double.isFinite(avg) ? avg : 0.0;
                        })
                        .toArray();

                for (double avg : shellAverages) {
                    line.append(",").append(avg);
                }

                write(line.toString());
            }
        };
        sim.addExporter(logger);
// Claude told me to add this to ensure that the executor is shut down when the simulation ends
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

        // sim.export();
        sim.preview();
    }
}


//import packages

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

//defining the space
public class BSimChemicalFieldTest {

    public static void main(String[] args) {
        final int bounds = 110;
        BSim sim = new BSim();
        sim.setDt(0.002); // difference in time
        sim.setTimeFormat("0.000");
        sim.setSolid(true, true, true); // always set solid to true otherwise get reflection
        sim.setBound(bounds, bounds, bounds); // size of simulation
        sim.setSimulationTime(120);

        final int boxcount = 32;// how many grid points - set up volume space -email will have name of limit?

//        parameters to play around with, main is threshold, diffusivity interacts with bounds and Dt
        final double decayRate = 0.07;
        final double diffusivity = 1380;
        final double threshold = 3e8;
        final double quantityadd = 3e10;
        final double firstwave = 180;
        final double secondwave = 1440;


        final BSimChemicalField field = new BSimChemicalField(sim, new int[]{boxcount, boxcount, boxcount}, diffusivity, decayRate);

//        placing the bacteria in the grid space
        final Vector<BSimBacterium> bacteria = new Vector<>();
        int totbac = 4500;
        while (bacteria.size() < totbac) {
            double u = Math.random();
            double v = Math.random();
            double theta = 2 * Math.PI * u;
            double phi = Math.acos(2 * v - 1);
            double x0 = Math.sin(phi) * Math.cos(theta);
            double y0 = Math.sin(phi) * Math.sin(theta);//this is just a sphere smth about a cude root
            double z0 = Math.cos(phi);
            double rand = Math.cbrt(Math.random());
            double x00 = rand * x0 * 0.5 * bounds;
            double y00 = rand * y0 * 0.5 * bounds;
            double z00 = rand * z0 * 0.5 * bounds;

            BSimBacterium p = new BSimBacterium(sim, new Vector3d(0.5 * bounds + x00, 0.5 * bounds + y00, 0.5 * bounds + z00)) {
                public boolean isFiring = false; //whether the wave has fired or not as you only want one fire
                public double fireStartTime = -1;
                public double fireDuration = 0.1;

                public boolean hasFiredFirstWave = false;

                public boolean isFiringSecondWave = false;
                public double fireStartTime2 = -1;
                public double fireDuration2 = 0.1;

                public boolean hasFiredSecondWave = false;


                @Override// this is the fire difffuse fire
                public void action() {
                    super.action();
                    double t = sim.getTime();

                    // Fire rates per time step  --- not used in this code
                    double fireRate = quantityadd / (fireDuration / sim.getDt());
                    double fireRate2 = quantityadd / (fireDuration2 / sim.getDt());

                    // Initial spike in center  ---- fires 7 times, changing this could change speed of the wave
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

                    // First wave firing logic (only once) --- if the concentration is higher than a threshold it will fire
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

                    // Second wave firing logic (only once)  ----- second wave same as first except it keeps going/ plateus
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
            if (!p.intersection(bacteria)) bacteria.add(p); // stops bacteria intersecting with each other
        }

        sim.setTicker(new BSimTicker() { // updates field when time ticks
            public void tick() {
                for (BSimBacterium b : bacteria) b.action();
                field.update();
            }
        });

        sim.setDrawer(new BSimP3DDrawer(sim, 800, 800) { // colours
            @Override
            public void scene(PGraphics3D p3d) {
                draw(field, new Color(0, 255, 255), (float) (255 / quantityadd));
                //for (BSimBacterium p : bacteria) draw(p, Color.RED);
            }
        });

        String resultsDir = BSimUtils.generateDirectoryPath("./results/");
        BSimLogger logger = new BSimLogger(sim, resultsDir + "SphereShellIntensity_centrifugal_2.csv") {
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
                for (int r_vox = 0; r_vox < maxShells; r_vox++) {
                    double r_min = r_vox * voxelSize_um;
                    double r_max = (r_vox + 1) * voxelSize_um;
                    double sum = 0.0;
                    int count = 0;
                    for (int x = 0; x < boxcount; x++) {
                        for (int y = 0; y < boxcount; y++) {
                            for (int z = 0; z < boxcount; z++) {
                                double x_um = (x + 0.5) * voxelSize_um;
                                double y_um = (y + 0.5) * voxelSize_um;
                                double z_um = (z + 0.5) * voxelSize_um;
                                double dx = x_um - centerX_um;
                                double dy = y_um - centerY_um;
                                double dz = z_um - centerZ_um;
                                double r = Math.sqrt(dx * dx + dy * dy + dz * dz);
                                if (r >= r_min && r < r_max) {
                                    double c = field.getConc(x, y, z);
                                    if (Double.isFinite(c) && c < 1e20) {
                                        sum += c;
                                        count++;
                                    }
                                }
                            }
                        }
                    }
                    double avg = (count > 0) ? sum / count : 0.0;
                    if (!Double.isFinite(avg)) avg = 0.0;
                    line.append(",").append(avg);
                }
                write(line.toString());
            }
        };
        sim.addExporter(logger);
        //sim.export();

        sim.preview();  // choose here whether to preview or export
    }
}
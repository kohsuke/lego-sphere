package org.kohsuke.lego;

import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import static java.lang.Math.*;
import static java.lang.Math.sin;
import static java.lang.Math.cos;
import static java.lang.System.out;
import java.util.ArrayList;
import java.util.List;

public class App {
    private static final class UVW {
        public final double u,v,w;

        private UVW(double u, double v, double w) {
            this.u = u;
            this.v = v;
            this.w = w;
        }
    }
    private static final class XYZ {
        public final double x,y,z;

        private XYZ(double x, double y, double z) {
            this.x = x;
            this.y = y;
            this.z = z;
        }

        public XYZ toPointReflection() {
            return new XYZ(-x,-y,-z);
        }

        /**
         * Distance from the origin.
         */
        public double r() {
            return sqrt(sq(x)+sq(y)+sq(z));
        }
    }

    private static final class Quadrant {
        public final int u, v;

        public Quadrant(int u, int v) {
            this.u = u;
            this.v = v;
        }

        public UVW transform(UVW p) {
            return new UVW(p.u*u,p.v*v,p.w);
        }
    }

    /**
     * Represents 4 quadrants.
     */
    private static Quadrant[] QUADRANTS = new Quadrant[] {
            new Quadrant( 1, 1),
            new Quadrant( 1,-1),
            new Quadrant(-1,1),
            new Quadrant(-1,-1)
    };

    /**
     * Longitude (theta) and latitude (phi)
     */
    private static final class PhiAndTheta {
        public final double phi,the;
        public final Color color;

        private PhiAndTheta(double phi, double the, Color color) {
            this.phi = phi;
            this.the = the;
            this.color = color;
        }
    }

    /**
     * A plate represented as a square whose four corners are in (theta,phi)
     */
    private static final class Plate {
        public final PhiAndTheta[] corners;

        private Plate(PhiAndTheta[] corners) {
            this.corners = corners;
        }
    }

    /**
     * Maps the piece-local (u,v,w) coordinate to the global (x,y,z) coordinate.
     */
    private abstract static class Piece {
        public final Color color;

        protected Piece(Color color) {
            this.color = color;
        }

        abstract XYZ map(UVW p);
    }

    /**
     * 3 different pieces. The other 3 are the point reflection.
     */
    private static Piece[] PIECES = new Piece[] {
        new Piece(Color.RED) {
            XYZ map(UVW p) {
                return new XYZ(p.w,-p.u,p.v);
            }
        },
        new Piece(Color.CYAN) {
            XYZ map(UVW p) {
                return new XYZ(-p.u,p.v,p.w);
            }
        },
        new Piece(Color.MAGENTA) {
            XYZ map(UVW p) {
                return new XYZ(p.v,p.w,-p.u);
            }
        }
    };

    /**
     * Represents 6 pieces that consistute the entire sphere.
     */

    public static void main(String[] args) throws IOException {
        final double l = 12;

        // square of the sphere radius
        final double r2 = sq(l/2)*3;

        // to decide in/out in the center of the plate
        // float d=0.5f;

        // to decide in/out at the far corner
        // float d=1f;

        final double d = 0.8;

        final double dx=1*d,dy=1*d,dz=0.4*d;

        // which voxels are occupied?
        boolean[] points = new boolean[(int)(l*l*l*3)];

        int mx=0,my=0,mz; // integer bound of regision

        int iz=0;
        for (double z=l/2; z<l*sqrt(3)/2; z+=0.4,iz++) {
            int iy=0;
            for (double y=0; y<l/2; y++,iy++) {
                int ix=0;
                for (double x=0; x<l*sqrt(2)/2; x++,ix++) {
                    boolean plate = sq(x+dx)+sq(y+dy)+sq(z+dz) < r2;
                    if(plate)   out.print('#');
                    else        out.print(' ');
                    points[c(ix,iy,iz,l)] = plate;
                }
                out.println();
                mx=ix;
            }
            out.println("-----------------");
            my=iy;
        }
        mz=iz;

        // compute the height map
        for (int z=0; z<mz; z++) {
            for (int y=0; y<my; y++) {
                for (int x=0; x<mx; x++) {
                    int h=0;
                    for (int i=z; i<mz; i++)
                        if(points[c(x,y,i,l)])
                            h++;
                    if(h==0)    out.print(' ');
                    else
                    if(h<6)     out.print((char)('0'+h));
                    else        out.print('#');
                }
                out.println();
            }
            out.println("-----------------");
        }

        List<PhiAndTheta> projections = new ArrayList<PhiAndTheta>();

        {// now let's compute the projection
            for (int v=0; v <my; v++) {
                for (int u=0; u<mx; u++) {
                    // compute the height at this (x,y)
                    int w=0;
                    for (int i=0; i<mz; i++)
                        if(points[c(u,v,i,l)])
                            w++;

                    if(w==0)    continue;

                    // center of this plate
                    final UVW p = new UVW(u+dx, v+dy, w*0.4+(l/2)+dz);

                    for (Quadrant q : QUADRANTS) {
//                        Piece piece = PIECES[1];
//                        int i=1;
//                        {{
                        for (Piece piece : PIECES) {
                            for(int i=0; i<2; i++) {
                                XYZ xyz = piece.map(q.transform(p));
                                if(i==1)
                                    xyz = xyz.toPointReflection();

                                if(xyz.y<0)
                                    continue; // this and adding 2 points to projection cancel out each other

                                // convert (x,y,z) to longitude and latitude
                                double phi   = asin(xyz.z/xyz.r());
                                double theta = acos(xyz.x/sqrt(sq(xyz.x)+sq(xyz.y)));

                                projections.add(new PhiAndTheta(phi,theta, piece.color));
                                projections.add(new PhiAndTheta(phi,theta+PI, piece.color));
                            }
                        }
                    }
                }
            }
        }

        System.out.println("projections="+projections.size());

        {// plot the computed projectoins on the miller projection map

            // dynamic range of y is (-yrange,yrange)
            double yrange = millerProjection(0.5*PI);

            BufferedImage img = new BufferedImage(2048,1502,BufferedImage.TYPE_INT_ARGB);
            Graphics2D g = img.createGraphics();
            for (PhiAndTheta p : projections) {
                int x = (int)(img.getWidth()*p.the/(2*PI));
                int y = (int)(img.getHeight()*(millerProjection(p.phi)+yrange)/(2*yrange));

                x=mod(x,img.getWidth());
                //y=mod(y,img.getHeight());


                g.setStroke(new BasicStroke(1));
                g.setColor(p.color);
                g.drawArc(x,y,10,10,0,360);
            }

            ImageIO.write(img,"PNG",new File("picture.png")); 
        }
    }

    private static double millerProjection(double phi) {
        return 1.25*sin(0.8*phi)/cos(0.8*phi);
    }

    private static int mod(int x, int mod) {
        return (x+mod)%mod;
    }

    private static int c(int x, int y, int z, double l) {
        int il = (int)l;
        return (z*il+y)*il+x;
    }

    /**
     * Computes x^2
     */
    private static double sq(double d) {
        return d*d;
    }

    /**
     * Computes sinh-1(x)
     */
    private static double sinhInv(double x) {
        return Math.log(x+sqrt(sq(x)+1));
    }
}

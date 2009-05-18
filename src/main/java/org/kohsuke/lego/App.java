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

        /**
         * Maps a piece-local (u,v,w) coordinate to the global (x,y,z) coordinate,
         * given the piece and the quadrant.
         */
        public XYZ toXYZ(Piece p) {
            return p.map(this);
        }

        /**
         * Computes the mirror point induced by the given quadrant.
         */
        public UVW flip(Quadrant q) {
            return q.transform(this);
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

        public PhiAndTheta toPhiAndTheta() {
            double phi   = asin(z/r());
            double theta = acos(x/sqrt(sq(x)+sq(y)));
            return new PhiAndTheta(phi,theta);
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

        private PhiAndTheta(double phi, double the) {
            this.phi = phi;
            this.the = the;
        }

        /**
         * Turn 180 degree along theta.
         */
        public PhiAndTheta rotatePi() {
            return new PhiAndTheta(phi,the+PI);
        }
    }

    /**
     * A plate represented as a square whose four corners are in (theta,phi)
     */
    private static final class Plate {
        public final PhiAndTheta[] corners;
        public final Color color;

        private Plate(Color color, PhiAndTheta... corners) {
            this.color =  color;
            this.corners = corners;
        }

        /**
         * Creates a plate that's turned 180 degree along theta
         */
        public Plate rotatePi() {
            PhiAndTheta[] x = new PhiAndTheta[4];
            for (int i=0; i<4; i++)
                x[i] = corners[i].rotatePi();
            return new Plate(color,x);
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
    private static Piece[] PIECES = buildPieces();

    private static Piece[] buildPieces() {
        // 3 primary pieces
        Piece p1 = new Piece(Color.RED) {
            XYZ map(UVW p) {
                return new XYZ(p.w, -p.u, p.v);
            }
        };
        Piece p2 = new Piece(Color.CYAN) {
            XYZ map(UVW p) {
                return new XYZ(-p.u, p.v, p.w);
            }
        };
        Piece p3 = new Piece(Color.MAGENTA) {
            XYZ map(UVW p) {
                return new XYZ(p.v, p.w, -p.u);
            }
        };

        // the other 3 are point reflection
        class MirrorPiece extends Piece {
            final Piece base;

            MirrorPiece(Piece base) {
                super(base.color);
                this.base = base;
            }

            XYZ map(UVW p) {
                return base.map(p).toPointReflection();
            }
        }


        return new Piece[] {p1,p2,p3,new MirrorPiece(p1),new MirrorPiece(p2),new MirrorPiece(p3)};
    }

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

        List<Plate> projections = new ArrayList<Plate>();

        {// now let's compute the projection
            for (int v=0; v <my; v++) {
                for (int u=0; u<mx; u++) {
                    // compute the height at this (x,y)
                    int w=0;
                    for (int i=0; i<mz; i++)
                        if(points[c(u,v,i,l)])
                            w++;
                    if(w==0)    continue;   // there's no plate here

                    for (Quadrant q : QUADRANTS) {// the same plate is copied to all 4 quadrants
                        for (Piece piece : PIECES) {// for all 6 pieces
                            // center of plate
                            final UVW p = new UVW(u+0.5, v+0.5, w*0.4+(l/2)+0.6).flip(q);
                            // corners of the plate
                            UVW[] corners = new UVW[4];
                            for (int j=0; j<4; j++)
                                corners[j] = new UVW(u+(j/2),v+((j==0||j==3)?0:1),p.w).flip(q);

                            XYZ xyz = p.toXYZ(piece);
                            if(xyz.y<0)
                                continue; // this and adding 2 points to projection cancel out each other

                            Plate plate = new Plate(piece.color,
                                    corners[0].toXYZ(piece).toPhiAndTheta(),
                                    corners[1].toXYZ(piece).toPhiAndTheta(),
                                    corners[2].toXYZ(piece).toPhiAndTheta(),
                                    corners[3].toXYZ(piece).toPhiAndTheta());

                            projections.add(plate);
                            projections.add(plate.rotatePi());
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
            for (Plate plate : projections) {
                // map 4 corners to screen x,y
                Point[] cp = new Point[4];
                for (int i = 0; i < cp.length; i++) {
                    PhiAndTheta p = plate.corners[i];
                    int x = (int)(img.getWidth()*p.the/(2*PI));
//                    x=mod(x,img.getWidth());
                    int y = (int)(img.getHeight()*(millerProjection(p.phi)+yrange)/(2*yrange));
                    cp[i] = new Point(x,y);
                }

                g.setStroke(new BasicStroke(1));
                g.setColor(plate.color);
                for (int i = 0; i < cp.length; i++) {
                    Point s = cp[i];
                    Point e = cp[(i+1)%4];
                    g.drawLine(s.x,s.y,e.x,e.y);
                }
            }

            ImageIO.write(img,"PNG",new File("picture2.png")); 
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

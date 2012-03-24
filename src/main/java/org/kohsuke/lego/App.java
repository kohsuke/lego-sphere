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
import java.util.Arrays;
import java.util.List;

public class App {
    /**
     * Earth is divided into 6 equal-size pieces like how valleyball is put together.
     * Looking from Sun to Earth, each piece would look something like:
     *
     * ^ V
     * |
     * |
     *
     *       -----------
     *      /           \
     *     /             \
     *     |      x      |
     *     \             /
     *      \           /
     *       -----------
     *                              ---> U
     *
     *
     * (U,V,W) is a coordinate local to a piece. Origin is depicted in "x" above.
     * studs would be facing upward, and so is the W-axis.
     *
     * W=0 is at the center of earth.
     * The bottom of this piece is a flat surface (which faces the inner cube.) That surface is W=l/2
     */
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

    /**
     * (X,Y,Z) is a coordinate global to Earth. Its origin is the center of earth.
     *
     * Looking from Sun to Earth, Z axis is the depth, X increases to the right, Y increases to the north pole.
     */
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
    private static List<Quadrant> QUADRANTS = Arrays.asList(
            new Quadrant( 1, 1),
            new Quadrant( 1,-1),
            new Quadrant(-1, 1),
            new Quadrant(-1,-1));

    /**
     * Polar coordinate. Longitude (theta) and latitude (phi)
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
    private static final class Surface {
        public final PhiAndTheta[] corners;
        public final Color color;

        private Surface(Color color, PhiAndTheta... corners) {
            this.color =  color;
            this.corners = corners;
        }

        private Surface(Piece piece, UVW c1, UVW c2, UVW c3, UVW c4) {
            this.color =  piece.color;
            this.corners = new PhiAndTheta[4];
            corners[0] = c1.toXYZ(piece).toPhiAndTheta();
            corners[1] = c2.toXYZ(piece).toPhiAndTheta();
            corners[2] = c3.toXYZ(piece).toPhiAndTheta();
            corners[3] = c4.toXYZ(piece).toPhiAndTheta();
        }

        /**
         * Creates a plate that's turned 180 degree along theta
         */
        public Surface rotatePi() {
            PhiAndTheta[] x = new PhiAndTheta[4];
            for (int i=0; i<4; i++)
                x[i] = corners[i].rotatePi();
            return new Surface(color,x);
        }
    }

    /**
     * Maps the piece-local (u,v,w) coordinate to the global (x,y,z) coordinate.
     */
    private abstract static class Piece {
        /**
         * Used to render the piece
         */
        public final Color color;

        protected Piece(Color color) {
            this.color = color;
        }

        abstract XYZ map(UVW p);
    }

    /**
     * 6 pieces that form the entire Earth.
     *
     * 3 are different, and the other 3 are the point reflections.
     */
    private static List<Piece> PIECES = buildPieces();

    private static List<Piece> buildPieces() {
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


        return Arrays.asList(p1,p2,p3,new MirrorPiece(p1),new MirrorPiece(p2),new MirrorPiece(p3));
    }
    
    private static class SurfaceList extends ArrayList<Surface> {
        public void put(Surface s) {
            add(s);
            add(s.rotatePi());
        }
    }

    /**
     * Width of 1x1x1 brick
     */
    private static final double BW = 1;

    /**
     * Height of 1x1x1 plate
     */
    private static final double BH = 0.4;

    // sphere contains a cube inside. This is the length of the cube edge
    private static final double l = 12;

    /**
     * Represents 6 pieces that constitute the entire sphere.
     */

    public static void main(String[] args) throws IOException {

        // from there, compute the sphere radius and its square
        final double r  = l/2*sqrt(3);
        final double r2 = sq(l/2)*3;

        // to decide in/out in the center of the plate
        final double d=0.8;
        // float double d=1;  // to decide in/out at the far corner
        // final double d = 0.8;

        final double dx=BW*d,dy=BW*d,dz=BH*d;

        // which voxels are occupied?
        boolean[] points = new boolean[(int)(l*l*l*3)];

        int mx=0,my=0,mz; // integer bound of region

        // looking from Sun into the center of earth (thus studs facing your way)
        // Z-axis is depth, X-axis is width, and Y-axis is height
        //
        // imagine the grid space of 1x1x1 plates. We scan them and decide which plate should be filled and which shouldn't
        int iz=0;
        for (double z=l/2; z<r; z+=BH,iz++) {
            int iy=0;
            for (double y=0; y<l/2; y+=BW,iy++) {
                int ix=0;
                for (double x=0; x<l*sqrt(2)/2; x+=BW,ix++) {
                    boolean plate = sq(x+dx)+sq(y+dy)+sq(z+dz) < r2;
                    if(plate)   out.print('#');
                    else        out.print(' ');
                    points[c(ix,iy,iz)] = plate;
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
                    // compute the number of 1x1x1 plate above z
                    int h=0;
                    for (int i=z; i<mz; i++)
                        if(points[c(x,y,i)])
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

        SurfaceList projections = new SurfaceList();

        {
            // now let's compute the projection
            // we scan a 4th of a piece, then flip them over to fill the entire piece, then copy them 6 times
            // to fully fill the entire earth.
            for (int v=0; v <my; v++) {
                for (int u=0; u<mx; u++) {
                    // compute the height at this (x,y)
                    int w=0;
                    for (int i=0; i<mz; i++)
                        if(points[c(u,v,i)])
                            w++;
                    if(w==0)    continue;   // there's no plate here

                    w--;    // we want w to be at the bottom of the plate, so that (u,v,w) points to the inner most point of the 1x1x1 plate

                    for (Quadrant q : QUADRANTS) {// the same plate is copied to all 4 quadrants
                        for (Piece piece : PIECES) {// for all 6 pieces
                            // 8 corners of the 1x1x1 plate at this position
                            UVW[][][] eight = eightCorners();
                            for (int _u=0; _u<2; _u++) {
                                for (int _v=0; _v<2; _v++) {
                                    for (int _w=0; _w<2; _w++) {
                                        eight[_u][_v][_w] = new UVW(u+_u,v+_v,(l/2)+(w+_w)*BH).flip(q);
                                    }
                                }                                
                            }

                            // 3 outer surfaces
                            projections.put(new Surface(piece,
                                    eight[0][0][1],
                                    eight[0][1][1],
                                    eight[1][1][1],
                                    eight[1][0][1]));

                            projections.put(new Surface(piece,
                                    eight[1][0][0],
                                    eight[1][0][1],
                                    eight[1][1][1],
                                    eight[1][1][0]));

                            projections.put(new Surface(piece,
                                    eight[0][1][0],
                                    eight[0][1][1],
                                    eight[1][1][1],
                                    eight[1][1][0]));
                        }
                    }
                }
            }
        }

        System.out.println("projections="+projections.size());

        {// plot the computed projections on the miller projection map

            // dynamic range of y is (-yrange,yrange)
            double yrange = millerProjection(0.5*PI);

            BufferedImage img = new BufferedImage(2048,1502,BufferedImage.TYPE_INT_ARGB);
            Graphics2D g = img.createGraphics();
            for (Surface plate : projections) {
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

    /**
     * Takes discrete integer (ix,iy,iz) position of a 1x1x1 plate cube, and converts that to the index in the array.
     */
    private static int c(int ix, int iy, int iz) {
        int il = (int)l;
        return (iz*il+iy)*il+ix;
    }

    /**
     * Allocates 2x2x2 array for holding 8 corners of 1x1x1 plate
     */
    private static UVW[][][] eightCorners() {
        UVW[][][] corners = new UVW[2][][];
        for (int i=0; i<2; i++) {
            corners[i] = new UVW[2][];
            for (int j=0; j<2; j++) {
                corners[i][j] = new UVW[2];
            }
        }

        return corners;
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

import edu.princeton.cs.algs4.*;

import java.awt.*;

public class SeamCarver {
    // create a seam carver object based on the given picture
    private Picture picture;

    public SeamCarver(Picture picture) {
        this.picture = picture;
    }

    // current picture
    public Picture picture() {
        return this.picture;
    }

    // width of current picture
    public int width() {
        return picture.width();
    }

    // height of current picture
    public int height() {
        return picture.height();
    }

    // energy of pixel at column x and row y
    public double energy(int x, int y) {
        if (x < 0 || x > this.width() || y < 0 || y > this.height())
            throw new IllegalArgumentException();

        if (x == 0 || y == 0 || x == this.width() - 1 || y == this.height() - 1)
            return 1000.0;

        double squareGradientX = calculateSquareGradientX(x, y);
        double squareGradientY = calculateSquareGradientY(x, y);

        return Math.sqrt(squareGradientX + squareGradientY);
    }

    private double calculateSquareGradientX(int x, int y) {
        Color colorPrev = new Color(this.picture.getRGB(x - 1, y));
        Color colorAfter = new Color(this.picture.getRGB(x + 1, y));

        double Rx = colorAfter.getRed() - colorPrev.getRed();
        double Gx = colorAfter.getGreen() - colorPrev.getGreen();
        double Bx = colorAfter.getBlue() - colorPrev.getBlue();

        return Rx * Rx + Gx * Gx + Bx * Bx;
    }

    private double calculateSquareGradientY(int x, int y) {
        Color colorPrev = new Color(this.picture.getRGB(x, y - 1));
        Color colorAfter = new Color(this.picture.getRGB(x, y + 1));

        double Ry = colorAfter.getRed() - colorPrev.getRed();
        double Gy = colorAfter.getGreen() - colorPrev.getGreen();
        double By = colorAfter.getBlue() - colorPrev.getBlue();

        return Ry * Ry + Gy * Gy + By * By;
    }

    // sequence of indices for horizontal seam
    public int[] findHorizontalSeam() {

        EdgeWeightedDigraph G = buildHorizontalGraph();

        AcyclicSP sp = new AcyclicSP(G, 0);
        int sourceMin = 0;
        int vMin = findMinPath(sp, sourceMin);
        int u;
        double distMin = sp.distTo(vMin);
        double dist;

        for (int i = 1; i < height(); i++) {
            sp = new AcyclicSP(G, from2Dto1D(i, 0));
            u = findMinPath(sp, from2Dto1D(i, 0));
            if (u == -1)
                continue;
            dist = sp.distTo(u);
            if (dist < distMin) {
                sourceMin = i;
                vMin = u;
                distMin = dist;
            }
        }

        sp = new AcyclicSP(G, sourceMin);


        int[] path = new int[width()];
        int i = 0;
        for (DirectedEdge edge : sp.pathTo(vMin)) {
            if (i == 0)
                path[i++] = edge.from();
            path[i++] = edge.to();
        }

        return path;
    }

    private int findMinPath(AcyclicSP sp, int source) {
        double min = sp.distTo(this.width() - 1);
        double end = width() - 1;
        double dist;
        int v = width() - 1;
        for (int i = 1; i < height(); i++) {
            dist = sp.distTo(from2Dto1D(i, width() - 1));
            if (dist < min)
                v = from2Dto1D(width() - 1, i);
        }
        return v;
    }

    private EdgeWeightedDigraph buildHorizontalGraph() {
        EdgeWeightedDigraph G = new EdgeWeightedDigraph(this.width() * this.height());

        for (int j = 0; j < this.width() - 1; j++) {
            for (int i = 0; i < this.height(); i++) {

                if (from2Dto1D(i - 1, j + 1) != -1)
                    G.addEdge(new DirectedEdge(
                            from2Dto1D(i, j),
                            from2Dto1D(i - 1, j + 1),
                            Math.abs(energy(j, i) + energy(j + 1, i - 1))
                    ));
                if (from2Dto1D(i + 1, j + 1) != -1)
                    G.addEdge(new DirectedEdge(
                            from2Dto1D(i, j),
                            from2Dto1D(i + 1, j + 1),
                            Math.abs(energy(j, i) + energy(j + 1, i + 1))
                    ));

                G.addEdge(new DirectedEdge(
                        from2Dto1D(i, j),
                        from2Dto1D(i, j + 1),
                        Math.abs(energy(j, i) + energy(j + 1, i))
                ));
            }
        }

        return G;
    }

    private int from2Dto1D(int x, int y) {
        if (x < 0 || x >= this.height() || y < 0 || y >= this.width())
            return -1;
        return x * this.width() + y;
    }

    // sequence of indices for vertical seam
    public int[] findVerticalSeam() {
        return new int[2];
    }

    // remove horizontal seam from current picture
    public void removeHorizontalSeam(int[] seam) {
        if (seam == null)
            throw new IllegalArgumentException();
    }

    // remove vertical seam from current picture
    public void removeVerticalSeam(int[] seam) {
        if (seam == null)
            throw new IllegalArgumentException();
        if (this.width() <= 1 || this.height() <= 1)
            throw new IllegalArgumentException();
    }

    public static void main(String[] args) {
        Picture picture1 = new Picture("6x5.png");
        SeamCarver sc = new SeamCarver(picture1);

        int[] path = sc.findHorizontalSeam();
        for (int i = 0; i < path.length; i++) {
            StdOut.println(path[i]);
        }
    }
}

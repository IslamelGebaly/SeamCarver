import edu.princeton.cs.algs4.Picture;

import java.awt.*;

public class SeamCarver {
    // create a seam carver object based on the given picture
    private Picture picture;
    private int width, height;

    public SeamCarver(Picture picture) {
        this.picture = picture;
        this.width = picture.width();
        this.height = picture.height();
    }

    // current picture
    public Picture picture() {
        return this.picture;
    }

    // width of current picture
    public int width() {
        return this.width;
    }

    // height of current picture
    public int height() {
        return this.height;
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
        return new int[2];
    }

    // sequence of indices for vertical seam
    public int[] findVerticalSeam() {
        double[] energy = new double[width() * height()];
        for (int j = 0; j < height(); j++) {
            for (int i = 0; i < width(); i++) {
                energy[from2Dto1D(i, j)] = energy(i, j);
            }
        }

        return new int[1];
    }

    private int from2Dto1D(int x, int y) {
        if (x < 0 || x >= this.height() || y < 0 || y >= this.width())
            return -1;
        return x * this.width() + y;
    }

    private int[][] adj(int i, int j) {
        int[][] adj;

        if (j > height() - 1)
            return null;
        if (i == 0 || i == width() - 1)
            adj = new int[2][2];
        else
            adj = new int[3][2];


        int n = i - 1;
        for (int index = 0; index < adj.length; ) {
            if (n != -1 && n != width()) {
                adj[index][0] = n;
                adj[index][1] = j + 1;
                index++;
            }
            n++;
        }

        return adj;
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

    }
}

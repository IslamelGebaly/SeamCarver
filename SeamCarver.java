import edu.princeton.cs.algs4.Picture;
import edu.princeton.cs.algs4.Stack;
import edu.princeton.cs.algs4.StdOut;

import java.awt.*;
import java.util.ArrayList;
import java.util.Arrays;

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

    private void setPicture(Picture picture) {
        this.picture = picture;
    }

    // width of current picture
    public int width() {
        return this.picture.width();
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
        Picture transpose = new Picture(height(), width());
        Picture temp = picture();
        for (int i = 0; i < picture.width(); i++) {
            for (int j = 0; j < picture.height(); j++)
                transpose.set(j, i, this.picture().get(i, j));
        }

        this.setPicture(transpose);
        int[] shortestPath = findVerticalSeam();
        this.setPicture(temp);

        for (int i = 0; i < width(); i++) {
            StdOut.println(energy(i, shortestPath[i]));
        }

        return shortestPath;
    }

    // sequence of indices for vertical seam
    public int[] findVerticalSeam() {
        final int WIDTH = width();
        final int HEIGHT = height();

        double[][] energy = new double[WIDTH][HEIGHT];

        for (int i = 0; i < WIDTH; i++) {
            for (int j = 0; j < HEIGHT; j++) {
                energy[i][j] = energy(i, j);
            }
        }

        double[][] distTo = calculateDistances(WIDTH, HEIGHT, energy);

        int[] shortestPath = new int[HEIGHT];
        int min;
        for (int j = 1; j < shortestPath.length; j++) {
            min = 0;
            for (int i = 0; i < WIDTH; i++) {
                min = distTo[i][j] < distTo[min][j] ? i : min;
            }
            shortestPath[j] = min;
        }

        shortestPath[0] = shortestPath[1];
        return shortestPath;
    }

    // remove horizontal seam from current picture
    public void removeHorizontalSeam(int[] seam) {
        if (seam == null)
            throw new IllegalArgumentException();
        if (seam.length != width())
            throw new IllegalArgumentException("String Length = " + String.valueOf(seam.length));
        if (width() <= 1)
            throw new IllegalArgumentException();

        Picture newPic = new Picture(width(), height() - 1);
        for (int i = 0; i < width(); i++) {
            if (seam[i] >= width() || seam[i] < 0)
                throw new IllegalArgumentException();
            if (i < width() - 1) {
                if (Math.abs(seam[i] - seam[i + 1]) > 1)
                    throw new IllegalArgumentException("Discrepancy = " +
                            String.valueOf(seam[i]) + "-" + String.valueOf(seam[i + 1]));
            }

            for (int j = 0, k = 0; j < height(); j++) {
                if (j == seam[i])
                    continue;
                newPic.set(i, k++, picture.get(i, j));
            }
        }

        setPicture(newPic);
    }

    // remove vertical seam from current picture
    public void removeVerticalSeam(int[] seam) {
        if (seam == null)
            throw new IllegalArgumentException();
        if (seam.length != width()) {
            throw new IllegalArgumentException("String Length = " + String.valueOf(seam.length));
        }
        if (height() <= 1)
            throw new IllegalArgumentException();

        Picture newPic = new Picture(width() - 1, height());

        for (int j = 0; j < height(); j++) {
            if (seam[j] >= width() || seam[j] < 0)
                throw new IllegalArgumentException();
            if (j < width() - 1) {
                if (Math.abs(seam[j] - seam[j + 1]) > 1)
                    throw new IllegalArgumentException("Discrepancy = " +
                            String.valueOf(seam[j]) + "-" + String.valueOf(seam[j + 1]));
            }

            for (int i = 0, l = 0; i < width(); i++) {
                if (i == seam[j])
                    continue;
                newPic.set(l++, j, picture.get(i, j));
            }
        }

        this.setPicture(newPic);
    }

    private int[][] adj(int x, int y) {
        int[][] adjacentPixels;
        int nexLevel = y + 1;
        if (y == height() - 1)
            return new int[0][0];
        if (x == 0 || x == width() - 1)
            adjacentPixels = new int[2][2];
        else
            adjacentPixels = new int[3][2];

        for (int index = 0, pixel = x - 1; index < adjacentPixels.length; ) {
            if (pixel < 0) {
                pixel++;
                continue;
            }

            adjacentPixels[index][0] = pixel++;
            adjacentPixels[index++][1] = nexLevel;
        }

        return adjacentPixels;
    }

    private double relax(double[][] energy, double pixel, double adj, int adjX, int adjY) {
        if (adj == 0)
            return pixel + energy[adjX][adjY];
        return Math.min(adj, pixel + energy[adjX][adjY]);
    }

    private double[][] calculateDistances(int w, int h, double[][] energy) {
        double[][] distTo = new double[w][h];

        Stack<ArrayList<Integer>> stack = new Stack<>();
        ArrayList<Integer> pixel;
        int adjX, adjY, pixelX, pixelY;

        for (int source = 0; source < w; source++) {
            stack.push(new ArrayList<>(Arrays.asList(source, 0)));
            while (!stack.isEmpty()) {
                pixel = stack.pop();
                pixelX = pixel.get(0);
                pixelY = pixel.get(1);
                for (int[] adjPixel : adj(pixelX, pixelY)) {
                    adjX = adjPixel[0];
                    adjY = adjPixel[1];

                    distTo[adjX][adjY] = relax(energy, distTo[pixelX][pixelY], distTo[adjX][adjY], adjX, adjY);

                    if (adjPixel[1] < h - 1)
                        stack.push(new ArrayList<>(Arrays.asList(adjPixel[0], adjPixel[1])));
                }
            }
        }

        return distTo;
    }

    public static void main(String[] args) {

    }
}

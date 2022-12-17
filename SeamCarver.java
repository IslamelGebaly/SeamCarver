import edu.princeton.cs.algs4.Picture;

import java.awt.Color;


public class SeamCarver {
    // create a seam carver object based on the given picture
    private Picture picture;
    private double[][] distTo;
    private double[][] energy;
    private Point[][] edgeTo;

    private class Point {
        int x, y;

        public Point(int x, int y) {
            this.x = x;
            this.y = y;
        }
    }

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
        return this.picture.height();
    }

    // energy of pixel at column x and row y
    public double energy(int x, int y) {
        if (x < 0 || x >= this.width() || y < 0 || y >= this.height())
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
        final int WIDTH = width();
        final int HEIGHT = height();

        Picture flippedPicture = new Picture(HEIGHT, WIDTH);
        Picture temp = picture();
        for (int i = 0; i < WIDTH; i++) {
            for (int j = 0; j < HEIGHT; j++)
                flippedPicture.set(j, i, this.picture().get(i, j));
        }

        this.setPicture(flippedPicture);
        int[] shortestPath = findVerticalSeam();
        this.setPicture(temp);

        return shortestPath;
    }

    // sequence of indices for vertical seam
    public int[] findVerticalSeam() {
        final int WIDTH = width();
        final int HEIGHT = height();

        this.energy = new double[WIDTH][HEIGHT];
        this.edgeTo = new Point[WIDTH][HEIGHT];
        this.distTo = new double[WIDTH][HEIGHT];

        for (int i = 0; i < WIDTH; i++) {
            for (int j = 0; j < HEIGHT; j++) {
                this.energy[i][j] = energy(i, j);
                if (j == 0)
                    this.distTo[i][j] = 0;
                else
                    this.distTo[i][j] = Double.POSITIVE_INFINITY;
                this.edgeTo[i][j] = new Point(-1, -1);
            }
        }

        calculateDistances(WIDTH, HEIGHT);

        int[] shortestPath = new int[HEIGHT];
        double minDist = Double.POSITIVE_INFINITY;
        int minPos = 0;

        for (int i = 0; i < WIDTH; i++) {
            if (minDist > distTo[i][HEIGHT - 1]) {
                minDist = distTo[i][HEIGHT - 1];
                minPos = i;
            }
        }

        Point vertix = new Point(minPos, HEIGHT - 1);

        for (int j = HEIGHT - 1; j >= 0; j--) {
            shortestPath[j] = vertix.x;
            if (j > 0)
                vertix = this.edgeTo[vertix.x][vertix.y];
        }

        return shortestPath;
    }

    // remove horizontal seam from current picture
    public void removeHorizontalSeam(int[] seam) {
        final int WIDTH = width();
        final int HEIGHT = height();

        if (seam == null)
            throw new IllegalArgumentException();
        if (seam.length != WIDTH)
            throw new IllegalArgumentException("String Length = " + String.valueOf(seam.length));
        if (HEIGHT <= 1)
            throw new IllegalArgumentException();

        Picture newPic = new Picture(WIDTH, HEIGHT - 1);
        for (int i = 0; i < WIDTH; i++) {
            if (seam[i] >= HEIGHT || seam[i] < 0)
                throw new IllegalArgumentException();
            if (i < WIDTH - 1) {
                if (Math.abs(seam[i] - seam[i + 1]) > 1)
                    throw new IllegalArgumentException("Discrepancy = " +
                            String.valueOf(seam[i]) + "-" + String.valueOf(seam[i + 1]));
            }

            for (int j = 0, k = 0; j < HEIGHT; j++) {
                if (j == seam[i])
                    continue;
                newPic.set(i, k++, picture.get(i, j));
            }
        }

        setPicture(newPic);
    }

    // remove vertical seam from current picture
    public void removeVerticalSeam(int[] seam) {
        final int WIDTH = width();
        final int HEIGHT = height();
        if (seam == null)
            throw new IllegalArgumentException();
        if (seam.length != HEIGHT) {
            throw new IllegalArgumentException("String Length = " + String.valueOf(seam.length));
        }
        if (WIDTH <= 1)
            throw new IllegalArgumentException();

        Picture newPic = new Picture(WIDTH - 1, HEIGHT);

        for (int j = 0; j < HEIGHT; j++) {
            if (seam[j] >= WIDTH || seam[j] < 0)
                throw new IllegalArgumentException();
            if (j < HEIGHT - 1) {
                if (Math.abs(seam[j] - seam[j + 1]) > 1)
                    throw new IllegalArgumentException("Discrepancy = " +
                            String.valueOf(seam[j]) + "-" + String.valueOf(seam[j + 1]));
            }

            for (int i = 0, l = 0; i < WIDTH; i++) {
                if (i == seam[j])
                    continue;
                newPic.set(l++, j, picture.get(i, j));
            }
        }

        this.setPicture(newPic);
    }

    private Point[] adj(Point pixel) {
        Point[] adj = {
                new Point(pixel.x - 1, pixel.y + 1),
                new Point(pixel.x, pixel.y + 1),
                new Point(pixel.x + 1, pixel.y + 1)
        };

        return adj;
    }

    private void relax(Point pixel, Point adj) {
        if (this.distTo[adj.x][adj.y] > this.distTo[pixel.x][pixel.y] + this.energy[adj.x][adj.y]) {
            this.edgeTo[adj.x][adj.y] = pixel;
            this.distTo[adj.x][adj.y] = this.distTo[pixel.x][pixel.y] + this.energy[adj.x][adj.y];
        }
    }

    private void calculateDistances(int w, int h) {
        Point pixel;
        Point[] adjacent;
        double minDistance = Double.POSITIVE_INFINITY;

        for (int j = 0; j < h - 1; j++) {
            for (int i = 0; i < w; i++) {
                pixel = new Point(i, j);
                adjacent = adj(pixel);
                if (adjacent[0].x > 0) {
                    relax(pixel, adjacent[0]);
                    minDistance = Math.min(this.distTo[adjacent[0].x][adjacent[0].y], minDistance);
                }
                relax(pixel, adjacent[1]);
                minDistance = Math.min(this.distTo[adjacent[1].x][adjacent[1].y], minDistance);
                if (adjacent[2].x < w) {
                    relax(pixel, adjacent[2]);
                    minDistance = Math.min(this.distTo[adjacent[2].x][adjacent[2].y], minDistance);
                }
            }
        }
    }

    public static void main(String[] args) {

    }
}

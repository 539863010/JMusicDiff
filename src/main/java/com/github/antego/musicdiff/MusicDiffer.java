package com.github.antego.musicdiff;

import javazoom.spi.mpeg.sampled.file.MpegAudioFileReader;

import javax.sound.sampled.*;
import java.io.*;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Objects;


public class MusicDiffer {
    private static int CHUNK_SIZE = (int)Math.pow(2, 10) * 4;
    private static final int[] RANGE = new int[] { 40, 80, 120, 180, 300 };
    private static final int FUZ_FACTOR = 2;


    public static void main(String[] args) throws Exception {
        byte[] audio = getRawAudioBytes(Paths.get(args[0]));
        int totalSize = audio.length;
        //todo round to bigger
        int chunksCount = totalSize / CHUNK_SIZE;
        Complex[][] result = new Complex[chunksCount][];

        for(int j = 0; j < chunksCount; j++) {
            Complex[] complexArray = new Complex[CHUNK_SIZE];

            for(int i = 0; i < CHUNK_SIZE; i++) {
                complexArray[i] = new Complex(audio[(j * CHUNK_SIZE) + i], 0);
            }

            result[j] = fft(complexArray);
        }

        getFootprint(result);

    }

    public static double[][] getFootprint(Complex[][] result) {
        double[][] highscores = new double[result.length][RANGE.length];
        for (int t = 0; t < result.length; t++) {
            for (int freq = RANGE[0]; freq < RANGE[RANGE.length - 1]; freq++) {
                double mag = Math.log(result[t][freq].abs() + 1);

                int index = getIndex(freq);

                if (mag > highscores[t][index]) {
                    highscores[t][index] = mag;
                }
            }
        }

        return highscores;
    }


    private static long hash(long p1, long p2, long p3, long p4) {
        return (p4 - (p4 % FUZ_FACTOR)) * 100000000
                + (p3 - (p3 % FUZ_FACTOR)) * 100000
                + (p2 - (p2 % FUZ_FACTOR)) * 100
                + (p1 - (p1 % FUZ_FACTOR));
    }

    public static int getIndex(int freq) {
        int i = 0;
        while (RANGE[i] < freq)
            i++;
        return i;
    }

    public static byte[] getRawAudioBytes(Path path) throws Exception {
        AudioInputStream in = new MpegAudioFileReader().getAudioInputStream(path.toFile());
        AudioFormat baseFormat = in.getFormat();
        AudioFormat decodedFormat =
                new AudioFormat(AudioFormat.Encoding.PCM_SIGNED,
                        baseFormat.getSampleRate(),
                        16,
                        baseFormat.getChannels(),
                        baseFormat.getChannels() * 2,
                        baseFormat.getSampleRate(),
                        false);
        AudioInputStream din = AudioSystem.getAudioInputStream(decodedFormat, in);

        byte[] data = new byte[4096];
        SourceDataLine line = getLine(decodedFormat);
        ByteArrayOutputStream out = new ByteArrayOutputStream();
        if (line != null) {
            int nBytesRead = 0;
            line.start();
            while (nBytesRead != -1) {
                nBytesRead = din.read(data, 0, data.length);
                if (nBytesRead > 0)
                    out.write(data, 0, nBytesRead);
            }
            line.stop();
            line.close();
            din.close();
        }
        in.close();
        return out.toByteArray();
    }

    private static SourceDataLine getLine(AudioFormat audioFormat)
            throws LineUnavailableException
    {
        SourceDataLine res = null;
        DataLine.Info info =
                new DataLine.Info(SourceDataLine.class, audioFormat);
        res = (SourceDataLine) AudioSystem.getLine(info);
        res.open(audioFormat);
        return res;
    }


    public static Complex[] fft(Complex[] x) {
        int n = x.length;

        // base case
        if (n == 1) return new Complex[] { x[0] };

        // radix 2 Cooley-Tukey FFT
        if (n % 2 != 0) { throw new RuntimeException("n is not a power of 2"); }

        // fft of even terms
        Complex[] even = new Complex[n/2];
        for (int k = 0; k < n/2; k++) {
            even[k] = x[2*k];
        }
        Complex[] q = fft(even);

        // fft of odd terms
        Complex[] odd  = even;  // reuse the array
        for (int k = 0; k < n/2; k++) {
            odd[k] = x[2*k + 1];
        }
        Complex[] r = fft(odd);

        // combine
        Complex[] y = new Complex[n];
        for (int k = 0; k < n/2; k++) {
            double kth = -2 * k * Math.PI / n;
            Complex wk = new Complex(Math.cos(kth), Math.sin(kth));
            y[k]       = q[k].plus(wk.times(r[k]));
            y[k + n/2] = q[k].minus(wk.times(r[k]));
        }
        return y;
    }


    public static class Complex {
        private final double re;   // the real part
        private final double im;   // the imaginary part

        // create a new object with the given real and imaginary parts
        public Complex(double real, double imag) {
            re = real;
            im = imag;
        }

        // return a string representation of the invoking Complex object
        public String toString() {
            if (im == 0) return re + "";
            if (re == 0) return im + "i";
            if (im <  0) return re + " - " + (-im) + "i";
            return re + " + " + im + "i";
        }

        // return abs/modulus/magnitude
        public double abs() {
            return Math.hypot(re, im);
        }

        // return angle/phase/argument, normalized to be between -pi and pi
        public double phase() {
            return Math.atan2(im, re);
        }

        // return a new Complex object whose value is (this + b)
        public Complex plus(Complex b) {
            Complex a = this;             // invoking object
            double real = a.re + b.re;
            double imag = a.im + b.im;
            return new Complex(real, imag);
        }

        // return a new Complex object whose value is (this - b)
        public Complex minus(Complex b) {
            Complex a = this;
            double real = a.re - b.re;
            double imag = a.im - b.im;
            return new Complex(real, imag);
        }

        // return a new Complex object whose value is (this * b)
        public Complex times(Complex b) {
            Complex a = this;
            double real = a.re * b.re - a.im * b.im;
            double imag = a.re * b.im + a.im * b.re;
            return new Complex(real, imag);
        }

        // return a new object whose value is (this * alpha)
        public Complex scale(double alpha) {
            return new Complex(alpha * re, alpha * im);
        }

        // return a new Complex object whose value is the conjugate of this
        public Complex conjugate() {
            return new Complex(re, -im);
        }

        // return a new Complex object whose value is the reciprocal of this
        public Complex reciprocal() {
            double scale = re*re + im*im;
            return new Complex(re / scale, -im / scale);
        }

        // return the real or imaginary part
        public double re() { return re; }
        public double im() { return im; }

        // return a / b
        public Complex divides(Complex b) {
            Complex a = this;
            return a.times(b.reciprocal());
        }

        // return a new Complex object whose value is the complex exponential of this
        public Complex exp() {
            return new Complex(Math.exp(re) * Math.cos(im), Math.exp(re) * Math.sin(im));
        }

        // return a new Complex object whose value is the complex sine of this
        public Complex sin() {
            return new Complex(Math.sin(re) * Math.cosh(im), Math.cos(re) * Math.sinh(im));
        }

        // return a new Complex object whose value is the complex cosine of this
        public Complex cos() {
            return new Complex(Math.cos(re) * Math.cosh(im), -Math.sin(re) * Math.sinh(im));
        }

        // return a new Complex object whose value is the complex tangent of this
        public Complex tan() {
            return sin().divides(cos());
        }



        // a static version of plus
        public static Complex plus(Complex a, Complex b) {
            double real = a.re + b.re;
            double imag = a.im + b.im;
            Complex sum = new Complex(real, imag);
            return sum;
        }

        // See Section 3.3.
        public boolean equals(Object x) {
            if (x == null) return false;
            if (this.getClass() != x.getClass()) return false;
            Complex that = (Complex) x;
            return (this.re == that.re) && (this.im == that.im);
        }

        // See Section 3.3.
        public int hashCode() {
            return Objects.hash(re, im);
        }
    }
}

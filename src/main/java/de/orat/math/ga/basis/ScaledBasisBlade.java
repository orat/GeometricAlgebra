// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


/**
 * ScaledBasisBlade.java
 *
 * Created on February 1, 2005, 11:41 AM
 *
 * Copyright 2005-2007 Daniel Fontijne, University of Amsterdam
 * fontijne@science.uva.nl
 *
 */
package de.orat.math.ga.basis;

import de.orat.math.ga.metric.MetricException;
import de.orat.math.ga.util.Bits;
import de.orat.math.ga.util.DoubleU;
import java.util.ArrayList;
import java.util.List;

/**
 * A simple class to represent a scaled basis blade.
 *
 * <p>mutable :( Should have made it immutable . . .
 *
 * <p>Could use subspace.util.Bits.lowestOneBit() and such to make
 * loops slightly more efficient, but didn't to keep text simple for the book.
 *
 * @author  fontijne
 */
/*final*/ public class ScaledBasisBlade implements Cloneable, InnerProductTypes {

    /** 
     * temp for testing.
     * 
     * @param args
     */
    public static void main(String args[]) {
        ScaledBasisBlade e1 = new ScaledBasisBlade(1);
        ScaledBasisBlade no = new ScaledBasisBlade(2);
        ScaledBasisBlade ni = new ScaledBasisBlade(4);

        double[][] m = new double[][]{
                {1.0, 0.0, 0.0, 0.0, 0.0},
                {0.0, 1.0, 0.0, 0.0, 0.0},
                {0.0, 0.0, 1.0, 0.0, 0.0},
                {0.0, 0.0, 0.0 ,0.0, -1.0},
                {0.0, 0.0, 0.0 , -1.0, 0.0}
        };

        try {
            Metric M = new Metric(m); {
            for (int i = 0; i < 32; i++)
                for (int j = 0; j < 32; j++) {
                    ScaledBasisBlade a = new ScaledBasisBlade(i);
                    ScaledBasisBlade b = new ScaledBasisBlade(j);
                    List<ScaledBasisBlade> R = Util.round(ScaledBasisBlade.gp(a, b, M), 1.0, 0.0001);
                }
            }
            /*
            ArrayList L = new ArrayList();
            long t = System.currentTimeMillis();
            for (int i = 0; i < 32; i++)
            for (int j = 0; j < 32; j++) {
                    ScaledBasisBlade a = new ScaledBasisBlade(i);
                    ScaledBasisBlade b = new ScaledBasisBlade(j);
                    ArrayList R = gp(a, b, M);
                    System.out.println(a + " " + b + " = " + R);
                    L.add(R);
            }

            System.out.println("time: " + (System.currentTimeMillis() - t));*/
        } catch (MetricException E) {
                System.out.println("Exception: " + E);
        }
    }

    /** 
     * Constructs an instance of BasisBlade.
     * 
     * @param b index (starts with 0)
     * @param s value weight
     */
    public ScaledBasisBlade(int b, double s) {
        bitmap = b;
        scale = s;
    }

    /** 
     * Constructs an instance of a unit BasisBlade.
     * 
     * @param b 
     */
    public ScaledBasisBlade(int b) {
        bitmap = b;
        scale = 1.0;
    }

    /** 
     * constructs an instance of a scalar BasisBlade.
     * 
     * @param s 
     */
    public ScaledBasisBlade(double s) {
        bitmap = 0;
        scale = s;
    }

    /** 
     * constructs an instance of a zero BasisBlade 
     */
    public ScaledBasisBlade() {
        bitmap = 0;
        scale = 0.0;
    }

    /**
     * Returns sign change due to putting the blade blades represented
     * by <code>a<code> and <code>b</code> into canonical order.
     * 
     * @param a
     * @param b
     * @return 
     */
    public static double canonicalReorderingSign(int a, int b) {
        // Count the number of basis vector flips required to
        // get a and b into canonical order.
        a >>= 1;
        int sum = 0;
        while (a != 0) {
            sum += Bits.bitCount(a & b);
            a >>= 1;
        }

        // even number of flips -> return 1
        // odd number of flips -> return 1
        return ((sum & 1) == 0) ? 1.0 : -1.0;
    }

    /** 
     * shortcut to outerProduct().
     * 
     * @param a
     * @param b
     * @return  
     */
    public static ScaledBasisBlade op(ScaledBasisBlade a, ScaledBasisBlade b) {
	return outerProduct(a, b);
    }

    /**
     * Outer/Wedge product.
     * 
     * @param a 
     * @param b
     * @return the outer product of two basis blades
     */
    public static ScaledBasisBlade outerProduct(ScaledBasisBlade a, ScaledBasisBlade b) {
	return gp_op(a, b, true);
    }

    /** 
     * shortcut to gp().
     * 
     * @param a
     * @param b
     * @return  
     */
    public static ScaledBasisBlade gp(ScaledBasisBlade a, ScaledBasisBlade b) {
	return geometricProduct(a, b);
    }

    /** 
     * Geometric product.
     * 
     * @param a
     * @param b
     * @return  the geometric product of two basis blades.
     */
    public static ScaledBasisBlade geometricProduct(ScaledBasisBlade a, ScaledBasisBlade b) {
	return gp_op(a, b, false);
    }

    /**
     * Geometric or outer product.
     * 
     * @param a
     * @param b
     * @param outer
     * @return the geometric product or the outer product
     * of 'a' and 'b'.
     */
    protected static ScaledBasisBlade gp_op(ScaledBasisBlade a, ScaledBasisBlade b, boolean outer) {
        // if outer product: check for independence
        if (outer && ((a.bitmap & b.bitmap) != 0))
            return new ScaledBasisBlade(0.0);

        // compute the bitmap:
        int bitmap = a.bitmap ^ b.bitmap;

        // compute the sign change due to reordering:
        double sign = canonicalReorderingSign(a.bitmap, b.bitmap);

        // return result:
        return new ScaledBasisBlade(bitmap, sign * a.scale * b.scale);
    }

    /** 
     * GeeometricProduct.
     * 
     * @param a
     * @param b
     * @param m
     * @return 
     */
    /*public static ScaledBasisBlade gp(ScaledBasisBlade a, ScaledBasisBlade b, double[] m) {
	return gp(a, b, m);
    }*/

    /**
     * Reverse.
     * 
     * A reverse operator is an involution operator, which means that twice application
     * reproduces the original object.
     * 
     * @return reverse of this basis blade (always a newly constructed blade)
     */
    public ScaledBasisBlade reverse() {
        // multiplier = (-1)^(x(x-1)/2)
        return new ScaledBasisBlade(bitmap, minusOnePow((grade() * (grade() - 1)) / 2) * scale);
    }

    /**
     * Grade inversion.
     * 
     * This is also called "involution" in Dorst2007.
     * 
     * @return grade inversion/involution of this basis blade (always a newly constructed blade)
     */
    public ScaledBasisBlade gradeInversion() {
        // multiplier is (-1)^x
        return new ScaledBasisBlade(bitmap, minusOnePow(grade()) * scale);
    }

    /**
     * Clifford_conjugate
     * 
     * The conjugate operation is an involution, which means that twice application 
     * results in the original object.
     * 
     * The conjugate is equivalent to the reverse if the blades does not multiplay to -1.
     * 
     * The conjugate of a basis blade is its inverse.
     * 
     * @return clifford conjugate of this basis blade (always a newly constructed blade)
     */
    public ScaledBasisBlade cliffordConjugate() {
        // multiplier is ((-1)^(x(x+1)/2)
        return new ScaledBasisBlade(bitmap, minusOnePow((grade() * (grade() + 1)) / 2) * scale);
    }

    /**
     * gp_restricted_NE_metric
     * 
     * Computes the geometric product of two basis blades in limited non-Euclidean metric.
     * 
     * @param a
     * @param b
     * @param m is an array of doubles giving the metric for each basis vector.
     * @return geometric product
     */
    public static ScaledBasisBlade gp(ScaledBasisBlade a, ScaledBasisBlade b, double[] m) {
        // compute the geometric product in Euclidean metric:
        ScaledBasisBlade result = geometricProduct(a, b);

        // compute the meet (bitmap of annihilated vectors):
        int bitmap = a.bitmap & b.bitmap;

        // change the scale according to the metric:
        int i = 0;
        while (bitmap != 0) {
            if ((bitmap & 1) != 0) result.scale *= m[i];
            i++;
            bitmap = bitmap >> 1;
        }

        return result;
    }

    public static List<ScaledBasisBlade> gp(ScaledBasisBlade a, ScaledBasisBlade b, Metric M) {
	return geometricProduct(a, b, M);
    }

    /**
     * Computes the geometric product of two basis blades in arbitary non-Euclidean metric.
     * 
     * @param a
     * @param b
     * @param M is an instance of Metric giving the metric (and precomputed eigen vectors).
     * @return an ArrayList, because the result does not have to be a single ScaledBasisBlade anymore . . .
     */
    public static List<ScaledBasisBlade> geometricProduct(ScaledBasisBlade a, ScaledBasisBlade b, Metric M) {
        // convert argument to eigenbasis
        List<ScaledBasisBlade> A = M.toEigenbasis(a);
        List<ScaledBasisBlade> B = M.toEigenbasis(b);

        List<ScaledBasisBlade> result = new ArrayList();
        for (int i = 0; i < A.size(); i++) {
            for (int j = 0; j < B.size(); j++) {
                result.add(gp(A.get(i), B.get(j), M.getEigenMetric()));
            }
        }

        return M.toMetricBasis(simplify(result));
    }

    /**
     * Inner product.
     * 
     * @param a
     * @param b
     * @param type gives the type of inner product:
     * LEFT_CONTRACTION,RIGHT_CONTRACTION, HESTENES_INNER_PRODUCT or 
     * MODIFIED_HESTENES_INNER_PRODUCT.
     * @return the geometric product of two basis blades
     */
    public static ScaledBasisBlade ip(ScaledBasisBlade a, ScaledBasisBlade b, int type) {
	return innerProductFilter(a.grade(), b.grade(), geometricProduct(a, b), type);
    }

    /**
     * Computes the inner product of two basis blades in limited non-Euclidean metric.
     * 
     * @param a
     * @param b
     * @param m is an array of doubles giving the metric for each basis vector.
     * @param type gives the type of inner product:
     * LEFT_CONTRACTION,RIGHT_CONTRACTION, HESTENES_INNER_PRODUCT or MODIFIED_HESTENES_INNER_PRODUCT.
     * @return 
     */
    public static ScaledBasisBlade ip(ScaledBasisBlade a, ScaledBasisBlade b, double[] m, int type) {
	return innerProductFilter(a.grade(), b.grade(), gp(a, b, m), type);
    }

    /**
     * Computes the inner product of two basis blades in arbitary non-Euclidean metric.
     * 
     * @param a left side blade
     * @param b right side blade
     * @param M is an instance of Metric giving the metric (and precomputed eigen vectors).
     * @param type gives the type of inner product:
     * LEFT_CONTRACTION,RIGHT_CONTRACTION, HESTENES_INNER_PRODUCT or MODIFIED_HESTENES_INNER_PRODUCT.
     * @return an ArrayList, because the result does not have to be a single ScaledBasisBlade anymore . . .
     * 
     * <p>Todo: no need to return an ArrayList here? Result is always a single blade?
     */
    public static List<ScaledBasisBlade> ip(ScaledBasisBlade a, ScaledBasisBlade b, Metric M, int type) {
	return innerProductFilter(a.grade(), b.grade(), geometricProduct(a, b, M), type);
    }

    /** 
     * shortcut to innerProduct(...).
     * 
     * @param a
     * @param b
     * @param type
     * @return  
     */
    /*public static ScaledBasisBlade ip(ScaledBasisBlade a, ScaledBasisBlade b, int type) {
	return innerProduct(a, b, type);
    }*/

    /** 
     * shortcut to innerProduct(...).
     * 
     * @param a
     * @param b
     * @param m
     * @param type
     * @return  
     */
    /*public static ScaledBasisBlade ip(ScaledBasisBlade a, ScaledBasisBlade b, double[] m, int type) {
	return innerProduct(a, b, m, type);
    }*/

    /** 
     * shortcut to innerProduct(...).
     * 
     * @param a
     * @param b
     * @param M
     * @param type
     * @return  
     */
    /*public static List<ScaledBasisBlade> ip(ScaledBasisBlade a, ScaledBasisBlade b, Metric M, int type) {
	return innerProduct(a, b, M, type);
    }*/

    /** 
     * Get grade.
     * 
     * @return the grade of this blade.
     */
    public int grade() {
	return Bits.bitCount(bitmap);
    }

    /**
     * Extract coordinates.
     * 
     * (C) Oliver Rettig 8.3.2022
     * 
     * @param n dimension
     * @return array of length of dimension with extracted coordinates, 0-values included.
     */
    /*public double[] extractCoordinate(int n){
        double[] result = new double[n];
        int h = Bits.highestOneBit(bitmap);
        if (h == -1) {
            result[0] = scale;
        } else {
            int l = Bits.lowestOneBit(bitmap);
            int index = 0; // index im result-array
            int offset = 0;
            int bitCount = 0;
            for (int bit=l;bit<h;bit++){
                if ((bitmap & 1) == 1){
                    bitCount++;
                    index+=bitCount/(n-1)+bitCount;
                }
            }
            result[index] = scale;
        }
        
        return result;
    }*/
    
    /**
     * Get index of the blade in the multivector, starting with 0.
     * 
     * @param n dimension of the vector space
     * @return index of the basis blade in the multivector
     */
    /*public int index(int n){
        int index = 0;
        int h = Bits.highestOneBit(bitmap);
        if (h >= 0) {
            int prev_i = 0;
            int g = 0; // Zahl der bisher durchlaufenen gesetzten bits
            for (int i=0;i<h;i++){
                if ((bitmap & 1) == 1){
                    int d = i-prev_i;
                    index+=(Util.binominal(n, g)+1)*d+i;
                    g++;
                    prev_i =i;
                }
            }
        }
        return index;
    }*/
    
    /**
     * @param i
     * @return pow(-1, i) 
     */
    public static int minusOnePow(int i) {
	return ((i & 1) == 0) ? 1 : -1;
    }

    @Override
    public boolean equals(Object O) {
        if (O instanceof ScaledBasisBlade) {
            ScaledBasisBlade B = (ScaledBasisBlade)O;
            return ((B.bitmap == bitmap) && (B.scale == scale));
        } else return false;
    }

    @Override
    public int hashCode() {
	return new Integer(bitmap).hashCode() ^ new Double(scale).hashCode();
    }

    @Override
    public String toString() {
		return toString(null);
    }

  
    /**
     * @param bvNames The names of the basis vector (e1, e2, e3) are used when
     * not available
     * 
     * @return Human readable String representation of the blade.
     */
    public String toString(String[] bvNames) {
        StringBuilder result = new StringBuilder();
        int i = 1;
        int b = bitmap;
        while (b != 0) {
            if ((b & 1) != 0) {
                if (result.length() > 0) result.append("^");
                if ((bvNames == null) || (i > bvNames.length) || (bvNames[i-1] == null))
                    result.append("e" + i);
                else result.append(bvNames[i-1]);
            }
            b >>= 1;
            i++;
        }
        return (result.length() == 0) ? Double.toString(scale) : scale + "*" + result.toString();
    }

    public ScaledBasisBlade copy() {
	return (ScaledBasisBlade) clone();
    }

    @Override
    public Object clone() {
	return new ScaledBasisBlade(bitmap, scale);
    }

    private static List<ScaledBasisBlade> innerProductFilter(int ga, int gb, List<ScaledBasisBlade> R, int type) {
        List<ScaledBasisBlade> result = new ArrayList();
        for (int i = 0; i < R.size(); i++) {
            ScaledBasisBlade B = innerProductFilter(ga, gb, (ScaledBasisBlade)R.get(i), type);
            if (B.scale != 0.0) result.add(B);
        }
        return result;
    }

    /**
     * Applies the rules to turn a geometric product into an inner product.
     * 
     * @param ga Grade of argument 'a'
     * @param gb Grade of argument 'b'
     * @param r the basis blade to be filter
     * @param type the type of inner product required:
     * LEFT_CONTRACTION,RIGHT_CONTRACTION, HESTENES_INNER_PRODUCT or MODIFIED_HESTENES_INNER_PRODUCT
     * @return Either a 0 basis blade, or 'r'
     */
    private static ScaledBasisBlade innerProductFilter(int ga, int gb, ScaledBasisBlade r, int type) {
        switch(type) {
            case LEFT_CONTRACTION:
                if ((ga > gb) || (r.grade() != (gb-ga))) {
                    return new ScaledBasisBlade();
                } else {
                    return r;
                }
            case RIGHT_CONTRACTION:
                if ((ga < gb) || (r.grade() != (ga-gb))){
                    return new ScaledBasisBlade();
                } else {
                    return r;
                }
            case HESTENES_INNER_PRODUCT:
                if ((ga == 0) || (gb == 0)) return new ScaledBasisBlade();
                // drop through to MODIFIED_HESTENES_INNER_PRODUCT
            // dot product?
            case MODIFIED_HESTENES_INNER_PRODUCT:
                if (Math.abs(ga - gb) == r.grade()) {
                    return r;
                } else {
                    return new ScaledBasisBlade();
                }
            default:
                return null;
        }
    }

    static private class simplifyComperator implements java.util.Comparator {
        @Override
        public int compare(Object o1, Object o2) {
            ScaledBasisBlade B1 = (ScaledBasisBlade)o1;
            ScaledBasisBlade B2 = (ScaledBasisBlade)o2;
            return ((B1.bitmap < B2.bitmap) ? -1 :
            ((B1.bitmap > B2.bitmap) ? 1 : Double.compare(B1.scale, B2.scale)));
        }
    }

    /** 
     * Simplifies a list (sum) of BasisBlades (destroys ordering of A!)
     * 
     * @param A
     * @return  list of simplified BasisBlades.
     */
    static protected List<ScaledBasisBlade> simplify(List<ScaledBasisBlade> A) {
        if (A.isEmpty()) return A;
        java.util.Collections.sort(A, new simplifyComperator());
        List<ScaledBasisBlade> result = new ArrayList();
        ScaledBasisBlade current = (A.get(0)).copy();
        for (int i = 1; i < A.size(); i++) {
            ScaledBasisBlade b = A.get(i);
            if (b.bitmap == current.bitmap) {
                current.scale += b.scale;
            } else {
                if (current.scale != 0.0) result.add(current);
                current = b.copy();
            }
        }
        if (current.scale != 0.0) result.add(current);
        return result;
    }

    /**
     * Rounds the scalar part of <code>this</code> to the nearest multiple X of <code>multipleOf</code>,
     * if |X - what| <= epsilon.
     * 
     * This is useful when eigenbasis is used to perform products in arbitrary
     * metric, which leads to small roundof errors.sYou don't want to keep these roundof errors if your
     * are computing a multiplication table.*
     * 
     * @param multipleOf
     * @param epsilon
     * @return 
     */
    public ScaledBasisBlade round(double multipleOf, double epsilon) {
        double a = DoubleU.round(scale, multipleOf, epsilon);
        if (a != scale)
            return new ScaledBasisBlade(bitmap, a);
        else return this;
    }

    /**
     * This bitmap specifies what basis vectors are
     * present in this basis blade
     */
    public int bitmap;

    /**
     * The scale/weight of the basis blade is represented by this double
     */
    public double scale;

} // end of class ScaledBasisBlade







































// leave this comment such that when view in HTML, the final part of the class will display properly

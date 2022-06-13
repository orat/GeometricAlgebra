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
 * Multivector.java
 *
 * Created on October 10, 2005, 7:37 PM
 *
 * Copyright 2005-2007 Daniel Fontijne, University of Amsterdam
 * fontijne@science.uva.nl
 *
 * TODO
 * - vielleicht besser in ein api package verschieben
 * - die spaceDim() Methode könnte fehlschlagen, wenn nicht alls BasisBlades in
 *   in einem Multivektor Werte != 0 haben. 
 *
 */
package de.orat.math.ga.basis;

import de.orat.math.ga.metric.MetricException;
import de.orat.math.ga.util.Bits;
import de.orat.math.ga.util.MathU;
import static de.orat.math.ga.util.MathU.sinh;
import java.util.*;

/**
 * This class implements a sample multivector class along with
 * some basic GA operations. Very low performance.
 *
 * <p>mutable :( Should have made it immutable . . .
 * @author  fontijne
 */
public class Multivector implements Cloneable, InnerProductTypes {
    
    public static void main(String[] args) {
	// setup conformal algebra:
	String[] bvNames = {"no", "e1", "e2", "e3", "ni"};
	double[][] m = new double[][]{
	    {0.0, 0.0, 0.0, 0.0, -1.0},
	    {0.0, 1.0, 0.0, 0.0, 0.0},
	    {0.0, 0.0, 1.0, 0.0, 0.0},
	    {0.0, 0.0, 0.0 ,1.0, 0.0},
	    {-1.0, 0.0, 0.0 , 0.0, 0.0}
	};

	Metric M = null;
	try {
	    M = new Metric(m);
	} catch (MetricException ex) {}

	Multivector no = createBasisVector(0);
	Multivector ni = createBasisVector(4);
	Multivector e1 = createBasisVector(1);
	Multivector e2 = createBasisVector(2);
	Multivector e3 = createBasisVector(3);

	Multivector A = e1.add(e2.op(e3).op(e1));
	System.out.println("A = " + A);
	System.out.println(new MultivectorType(A));

	/*
	Multivector a = e1.add(e2.op(e3));
	Multivector ai1 = a.generalInverse(M);
	Multivector T1 = ai1.gp(a, M).sub(1.0).compress(1e-7);
	Multivector T2 = a.gp(ai1, M).sub(1.0).compress(1e-7);
	if (!T1.isNull()) {
	    System.out.println("More Wha 1 !!!");
	}
	if (!T2.isNull()) {
	    System.out.println("More Wha 2!!!");
	}*/


	/*int dim = 5;
	for (int i = 0; i < 1000; i++) {
	    if ((i % 10) == 0) System.out.println("i = " + i);
	    Multivector a = Multivector.getRandomVersor(dim, (int)(Math.random() * dim), 1.0);
	    Multivector ai1 = null;
	    Multivector ai2 = null;
	    try {
		ai1 = a.generalInverse(M);
	    } catch (Exception ex) {
	    }
	    try {
		ai2 = a.versorInverse(M);
	    } catch (Exception ex) {
	    }
	    if ((ai1 == null) ^ (ai2 == null)) {
		System.out.println("Wha!!!");
	    }
	    if (ai1 != null) {
		Multivector T1 = ai1.gp(a, M).sub(1.0).compress(1e-7);
		Multivector T2 = a.gp(ai1, M).sub(1.0).compress(1e-7);
		if (!T1.isNull()) {
		    Multivector X = a.gp(ai1, M).sub(1.0).compress(1e-7);
		    System.out.println("More Wha 1 !!!");
		}
		if (!T2.isNull()) {
		    System.out.println("More Wha 2!!!");
		}
//		System.out.println("T1: " + T1.toString(bvNames));
//		System.out.println("T2: " + T2.toString(bvNames));

	    }

	}*/

	/*
	Multivector B = new Multivector(-1.45);
	Multivector R = B.cos(M);
	Multivector R2 = B.cosSeries(M, 24);

	System.out.println("B = " + B.toString(bvNames) + ",");
	System.out.println("R1 = " + R.toString(bvNames) + ",");
	System.out.println("R2 = " + R2.toString(bvNames) + ",");

	B = e1.op(e3).gp(1.33334);
	R = B.cos(M);
	R2 = B.cosSeries(M, 24);

	System.out.println("B = " + B.toString(bvNames) + ",");
	System.out.println("R1 = " + R.toString(bvNames) + ",");
	System.out.println("R2 = " + R2.toString(bvNames) + ",");
	 */

    }


    /** 
     * Creates a new instance of Multivector.
     */
    public Multivector() {
	blades = new ArrayList();
    }

    /** 
     * Creates a new instance of Multivector.
     * 
     * @param s scalar value
     */
    public Multivector(double s) {
	this(new ScaledBasisBlade(s));
    }

    /** 
     * Creates a new instance of Multivector.
     * 
     * do not modify 'B' for it is not copied.
     * 
     * @param B list of scaled blades
     */
    public Multivector(List<ScaledBasisBlade> B) {
	blades = B;
    }

    /** 
     * Creates a new instance of Multivector.
     * 
     * do not modify 'B' for it is not copied.
     * 
     * @param B 
     */
    public Multivector(ScaledBasisBlade B) {
	blades = new ArrayList<>(); //ArrayListU.newArrayList(B);
        blades.add(B);
    }

    public Multivector copy() throws CloneNotSupportedException {
	return (Multivector) clone();
    }

    @Override
    public Object clone() throws CloneNotSupportedException {
        List<ScaledBasisBlade> clonedBasisBlades = new ArrayList<>(blades);
        Multivector M = new Multivector(/*(ArrayList*//*(List<BasisBlade>)gblades.clone()*/clonedBasisBlades);
        M.bladesSorted = this.bladesSorted;
        return M;
    }

    @Override
    public boolean equals(Object O) {
        if (O instanceof Multivector B) {
            Multivector zero = Multivector.this.sub(B);
            return (zero.blades.isEmpty());
        } else {
            return false;
        }
    }

    @Override
    public int hashCode() {
        int hc = 0;
        for (int i = 0; i < blades.size(); i++)
            hc ^= blades.get(i).hashCode();
        return hc;
    }

    @Override
    public String toString() {
	return toString(null);
    }

    /**
     * toString.
     * 
     * @param bvNames The names of the basis vector (e1, e2, e3) are used when
     * not available
     * @return String representation of the multivector
     */
    public String toString(String[] bvNames) {
        if (blades.isEmpty()) {
            return "0";
        } else {
            StringBuilder result = new StringBuilder();
            for (int i = 0; i < blades.size(); i++) {
                ScaledBasisBlade b = blades.get(i);
                String S = b.toString(bvNames);
                if (i == 0) {
                    result.append(S);
                } else if (S.charAt(0) == '-') {
                    result.append(" - ");
                    result.append(S.substring(1));
                } else {
                    result.append(" + ");
                    result.append(S);
                }
            }
            return result.toString();
        }
    }

    /**
     * Create a basis vector.
     * 
     * @param idx 
     * @return basis vector 'idx' range [0 ... dim)
     */
    protected static Multivector createBasisVector(int idx) {
	return new Multivector(new ScaledBasisBlade(1 << idx));
    }
    /**
     * Create a scaled basis vector.
     * 
     * @param idx
     * @param s scale
     * @return basis vector 'idx' range [0 ... dim)
     */
    public static Multivector createBasisVector(int idx, double s) {
        ScaledBasisBlade basisBlade = new ScaledBasisBlade(1 << idx,s);
	return new Multivector(basisBlade);
    }
    
    /**
     * Get random vector.
     * 
     * @param dim 
     * @param scale
     * @return 'dim'-dimensional random vector with coordinates in range [-scale, scale] 
     */
    public static Multivector getRandomVector(int dim, double scale) {
        ArrayList result = new ArrayList(dim);
        for (int i = 0; i < dim; i++)
            result.add(new ScaledBasisBlade(1 << i, 2 * scale * (Math.random() - 0.5)));
        return new Multivector(result);
    }

    /**
     * Get random blade.
     * 
     * @param dim
     * @param grade
     * @param scale
     * @return 'dim'-dimensional random blade with coordinates in range [-scale, scale]
     */
    public static Multivector getRandomBlade(int dim, int grade, double scale) {
        Multivector result = new Multivector(2 * scale * (Math.random() - 0.5));
        for (int i = 1; i <= grade; i++)
            result = result.op(getRandomVector(dim, scale));
        return result;
    }

    /**
     * Get random versor.
     * 
     * @param dim
     * @param grade
     * @param scale
     * @return 'dim'-dimensional random blade with coordinates in range [-scale, scale]
     */
    public static Multivector getRandomVersor(int dim, int grade, double scale) {
	return getRandomVersor(dim, grade, scale, null);
    }

    /**
     * Get random versor.
     * 
     * @param dim
     * @param grade
     * @param scale
     * @return 'dim'-dimensional random blade with coordinates in range [-scale, scale]
     * @param metric can either be null, Metric or double[]
     */
    public static Multivector getRandomVersor(int dim, int grade, double scale, Object metric) {
        Multivector result = new Multivector(2 * scale * (Math.random() - 0.5));
        for (int i = 1; i <= grade; i++)
                result = result.gp(getRandomVector(dim, scale), metric);
        return result;
    }

    /**
     * Geometric product with a scalar.
     * 
     * @param a scalar
     * @return geometric product of this with a scalar
     */
    public Multivector gp(double a) {
        if (a == 0.0) {
            return new Multivector();
        } else {
            ArrayList result = new ArrayList(blades.size());
            for (int i = 0; i < blades.size(); i++) {
                ScaledBasisBlade b = blades.get(i);
                result.add(new ScaledBasisBlade(b.bitmap, b.scale * a));
            }
            return new Multivector(result);
        }
    }

    /**
     * Geometric product.
     * 
     * @param x multivector
     * @return  geometric product of this with a 'x' without using the metric.
     */
    protected Multivector gp(Multivector x) {
        List<ScaledBasisBlade> result = new ArrayList(blades.size() * x.blades.size());

        // loop over basis blade of 'this'
        for (int i = 0; i < blades.size(); i++) {
            ScaledBasisBlade B1 = blades.get(i);

            // loop over basis blade of 'x'
            for (int j = 0; j < x.blades.size(); j++) {
                ScaledBasisBlade B2 = x.blades.get(j);
                result.add(ScaledBasisBlade.gp(B1, B2));
            }
        }
        return new Multivector(simplify(result));
    }

    /**
     * Geometric product.
     * 
     * @param x multivector
     * @param M metric
     * @return  geometric product of this with a 'x' using metric 'M'
     */
    public Multivector gp(Multivector x, Metric M) {
        List<ScaledBasisBlade> result = new ArrayList(blades.size() * x.blades.size());
        for (int i = 0; i < blades.size(); i++) {
            ScaledBasisBlade B1 = blades.get(i);
            for (int j = 0; j < x.blades.size(); j++) {
                ScaledBasisBlade B2 = x.blades.get(j);
                result.addAll(ScaledBasisBlade.gp(B1, B2, M));
            }
        }
        return new Multivector(simplify(result));
    }

    /**
     * Geometric product.
     * 
     * @param x 
     * @param m 
     * @return geometric product of this with a 'x' using metric 'm' 
     */
    protected Multivector gp(Multivector x, double[] m) {
        ArrayList result = new ArrayList(blades.size() * x.blades.size());
        for (int i = 0; i < blades.size(); i++) {
            ScaledBasisBlade B1 = blades.get(i);
            for (int j = 0; j < x.blades.size(); j++) {
                ScaledBasisBlade B2 = x.blades.get(j);
                result.add(ScaledBasisBlade.gp(B1, B2, m));
            }
        }
        return new Multivector(simplify(result));
    }

    /**
     * Outer product.
     * 
     * @param x 
     * @return outer product of this with 'x' 
     */
    protected Multivector op(Multivector x) {
        ArrayList result = new ArrayList(blades.size() * x.blades.size());
        for (int i = 0; i < blades.size(); i++) {
            ScaledBasisBlade B1 = blades.get(i);
            for (int j = 0; j < x.blades.size(); j++) {
                ScaledBasisBlade B2 = x.blades.get(j);
                result.add(ScaledBasisBlade.op(B1, B2));
            }
        }
        return new Multivector(simplify(result));
    }

    /**
     * Inner product.
     * 
     * @param x 
     * @param type gives the type of inner product:
     * LEFT_CONTRACTION,
     * RIGHT_CONTRACTION,
     * HESTENES_INNER_PRODUCT or
     * MODIFIED_HESTENES_INNER_PRODUCT.
     * @return inner product of this with a 'x'
     */
    protected Multivector ip(Multivector x, int type) {
        ArrayList<ScaledBasisBlade> result = new ArrayList(blades.size() * x.blades.size());
        for (int i = 0; i < blades.size(); i++) {
            ScaledBasisBlade B1 = blades.get(i);
            for (int j = 0; j < x.blades.size(); j++) {
                ScaledBasisBlade B2 = x.blades.get(j);
                result.add(ScaledBasisBlade.ip(B1, B2, type));
            }
        }
        return new Multivector(simplify(result));
    }

    /**
     * Inner product.
     * 
     * @param x
     * @param M
     * @param type gives the type of inner product:
     * LEFT_CONTRACTION,
     * RIGHT_CONTRACTION,
     * HESTENES_INNER_PRODUCT or
     * MODIFIED_HESTENES_INNER_PRODUCT.
     * @return inner product of this with a 'x' using metric 'M'
     */
    public Multivector ip(Multivector x, Metric M, int type) {
        List<ScaledBasisBlade> result = new ArrayList(blades.size() * x.blades.size());
        for (int i = 0; i < blades.size(); i++) {
            ScaledBasisBlade B1 = blades.get(i);
            for (int j = 0; j < x.blades.size(); j++) {
                ScaledBasisBlade B2 = x.blades.get(j);
                result.addAll(ScaledBasisBlade.ip(B1, B2, M, type));
            }
        }
        return new Multivector(simplify(result));
    }

    /**
     * Inner product.
     * 
     * @param x 
     * @param m
     * @param type gives the type of inner product:
     * LEFT_CONTRACTION,
     * RIGHT_CONTRACTION,
     * HESTENES_INNER_PRODUCT or
     * MODIFIED_HESTENES_INNER_PRODUCT.
     * @return inner product of this with a 'x' using metric 'm'
     */
    protected Multivector ip(Multivector x, double[] m, int type) {
        List<ScaledBasisBlade> result = new ArrayList(blades.size() * x.blades.size());
        for (int i = 0; i < blades.size(); i++) {
            ScaledBasisBlade B1 = blades.get(i);
            for (int j = 0; j < x.blades.size(); j++) {
                ScaledBasisBlade B2 = x.blades.get(j);
                result.add(ScaledBasisBlade.ip(B1, B2, m, type));
            }
        }
        return new Multivector(simplify(result));
    }

    /**
     * Scalar product.
     * 
     * @param x 
     * @return scalar product of this with a 'x' but without respecting a metric.
     */
    protected double scp(Multivector x) {
	return ip(x, LEFT_CONTRACTION).scalarPart();
    }

    /**
     * Scalar product.
     * 
     * @param x 
     * @param m
     * @return scalar product of this with a 'x' using metric 'm'
     */
    /*public double scalarProduct(Multivector x, double[] m) {
	return ip(x, m, LEFT_CONTRACTION).scalarPart();
    }*/

    /**
     * Scalar product.
     * 
     * @param x
     * @param M metric
     * @return scalar product of this with a 'x' using metric 'M'
     */
    public double scp(Multivector x, Metric M) {
	return ip(x, M, LEFT_CONTRACTION).scalarPart();
    }

    /** 
     * Scalar product.
     * 
     * @param x
     * @param m
     * @return  
     */
    protected double scp(Multivector x, double[] m) {
        return ip(x, m, LEFT_CONTRACTION).scalarPart();
	//return scalarProduct(x, m);
    }

    /**
     * Add a scalar.
     * 
     * @param a 
     * @return sum of this with scalar 'a' 
     */
    public Multivector add(double a) {
        List<ScaledBasisBlade> result = new ArrayList(blades.size() + 1);
        result.addAll(blades);
        result.add(new ScaledBasisBlade(a));
        return new Multivector(simplify(cloneArrayElements(result)));
    }

    /**
     * Addition.
     * 
     * @param x 
     * @return sum of this with 'x'
     */
    public Multivector add(Multivector x) {
        List<ScaledBasisBlade> result = new ArrayList(blades.size() + x.blades.size());
        result.addAll(blades);
        result.addAll(x.blades);
        return new Multivector(simplify(cloneArrayElements(result)));
    }

    /**
     * Subtraction.
     * 
     * @param a 
     * @return turn this - scalar 'a' 
     */
    public Multivector sub(double a) {
	return add(-a);
    }

    /**
     * Subtraction.
     * 
     * @param x 
     * @return this - 'x' 
     */
    public Multivector sub(Multivector x) {
        List<ScaledBasisBlade> result = new ArrayList(blades.size() + x.blades.size());
        result.addAll(blades);
        result.addAll(x.gp(-1.0).blades);
        return new Multivector(simplify(cloneArrayElements(result)));
    }

    /** 
     * Exponential.
     * 
     * @return exponential of this 
     */
    protected Multivector exp() {
	return exp(null, 12);
    }
    /**
     * Exponential.
     * 
     * @param M metric
     * @return exponential of this under metric 'M'
     */
    public Multivector exp(Metric M) {
	return exp(M, 12);
    }
    /**
     * @param m metric
     * @return exponential of this under metric 'm'
     */
    protected Multivector exp(double[] m) {
	return exp(m, 12);
    }
    /** 
     * evaluates exp(this) using special cases if possible, using series otherwise.
     * 
     * @param M
     * @param order
     * @return  
     */
    protected Multivector exp(Object M, int order) {
        // check out this^2 for special cases
        Multivector A2 = this.gp(this, M).compress();
        if (A2.isNull(1e-8)) {
            // special case A^2 = 0
            return this.add(1);
        } else if (A2.isScalar()) {
            double a2 = A2.scalarPart();
            // special case A^2 = +-alpha^2
            if (a2 < 0) {
                double alpha = Math.sqrt(-a2);
                return gp(Math.sin(alpha) / alpha).add(Math.cos(alpha));
            }
            //hey: todo what if a2 == 0?
            else {
                double alpha = Math.sqrt(a2);
                return gp(MathU.sinh(alpha) / alpha).add(MathU.cosh(alpha));
            }
        } else return expSeries(M, order);
    }

    /** 
     * Evaluates exp using series ...(== SLOW & INPRECISE!)
     * 
     * @param M
     * @param order
     * @return  
     */
    protected Multivector expSeries(Object M, int order) {
        // first scale by power of 2 so that its norm is ~ 1
        long scale=1; {
            double max = this.norm_e();
            if (max > 1.0) scale <<= 1;
            while (max > 1.0) {
                max = max / 2;
                scale <<= 1;
            }
        }

        Multivector scaled = this.gp(1.0 / scale);

        // taylor approximation
        Multivector result = new Multivector(1.0); {
            Multivector tmp = new Multivector(1.0);
            for (int i = 1; i < order; i++) {
                tmp = tmp.gp(scaled.gp(1.0 / i), M);
                result = result.add(tmp);
            }
        }

        // undo scaling
        while (scale > 1) {
            result = result.gp(result, M);
            scale >>>= 1;
        }
        return result;
    }


    /** 
     * @return sin of this 
     */
    protected Multivector sin() {
	return sin(null, 12);
    }
    /**
     * @param M 
     * @return sin of this under metric 'M'  
     */
    public Multivector sin(Metric M) {
        return sin(M, 12);
    }
    /**
     * @param m 
     * @return sin of this under metric 'm' 
     */
    protected Multivector sin(double[] m) {
	return sin(m, 12);
    }

    protected Multivector sin(Object M, int order) {
        // check out this^2 for special cases
        Multivector A2 = this.gp(this, M).compress();
        if (A2.isNull(1e-8)) {
            // special case A^2 = 0
            return this;
        } else if (A2.isScalar()) {
            double a2 = A2.scalarPart();
            // special case A^2 = +-alpha^2
            if (a2 < 0) {
                double alpha = Math.sqrt(-a2);
                return gp(sinh(alpha) / alpha);
            }
            //hey: todo what if a2 == 0?
            else {
                double alpha = Math.sqrt(a2);
                return gp(Math.sin(alpha) / alpha);
            }
        } else return sinSeries(M, order);
    }

    /** 
     * Evaluates sin using series ...(== SLOW & INPRECISE!)
     * 
     * @param M
     * @param order
     * @return  
     */
    protected Multivector sinSeries(Object M, int order) {
        Multivector scaled = this;
        // taylor approximation
        Multivector result = scaled;
        Multivector tmp = scaled;
        int sign = -1;
        for (int i = 2; i < order; i ++) {
            tmp = tmp.gp(scaled.gp(1.0 / i), M);
            if ((i & 1) != 0) {// only the odd part of the series
                result = result.add(tmp.gp((double)sign));
                sign *= -1;
            }
        }
        return result;
    }

    /** 
     * @return cos of this 
     */
    protected Multivector cos() {
        return cos(null, 12);
    }
    /**
     * @param M 
     * @return cos of this under metric 'M'
     */
    public Multivector cos(Metric M) {
	return cos(M, 12);
    }
    /**
     * @param m 
     * @return cos of this under metric 'm'
     */
    protected Multivector cos(double[] m) {
	return cos(m, 12);
    }

    protected Multivector cos(Object M, int order) {
        // check out this^2 for special cases
        Multivector A2 = this.gp(this, M).compress();
        if (A2.isNull(1e-8)) {
            // special case A^2 = 0
            return new Multivector(1);
        } else if (A2.isScalar()) {
            double a2 = A2.scalarPart();
            // special case A^2 = +-alpha^2
            if (a2 < 0) {
                    double alpha = Math.sqrt(-a2);
                    return new Multivector(MathU.cosh(alpha));
            }
            //hey: todo what if a2 == 0?
            else {
                    double alpha = Math.sqrt(a2);
                    return new Multivector(Math.cos(alpha));
            }
        } else return cosSeries(M, order);
    }

    /**
     * Evaluates cos using series ...(== SLOW & INPRECISE!).
     * 
     * @param M
     * @param order
     * @return  
     */
    protected Multivector cosSeries(Object M, int order) {
        Multivector scaled = this;
        // taylor approximation
        Multivector result = new Multivector(1.0); 
        Multivector tmp = scaled;
        int sign = -1;
        for (int i = 2; i < order; i ++) {
            tmp = tmp.gp(scaled.gp(1.0 / i), M);
            if ((i & 1) == 0) {// only the even part of the series
                result = result.add(tmp.gp((double)sign));
                sign *= -1;
            }
        }
        return result;
    }

    /**
     * Unit multivector under euclidian norm.
     * 
     * @return unit under Euclidean norm
     * @throws java.lang.ArithmeticException if multivector is null.
     */
    protected Multivector unit() {
	return unit_r();
    }
    /*public Multivector unit_e(Metric m) {
	return unit_r(m);
    }*/
    /**
     * Unit multivector under reverse norm.
     * 
     * @throws java.lang.ArithmeticException if multivector is null
     * @return unit under 'reverse' norm (this / sqrt(abs(this.reverse(this))))
     */
    protected Multivector unit_r() {
        double s = scp(reverse());
        if (s == 0.0) throw new java.lang.ArithmeticException("null multivector");
        else return this.gp(1d / Math.sqrt(Math.abs(s)));
    }
    /**
     * @throws java.lang.ArithmeticException if multivector is null
     * @param m
     * @return unit under 'reverse' norm (this / sqrt(abs(this.reverse(this))))
     */
    /*public Multivector unit_r(double[] m) {
        double s = scp(reverse(), m);
        if (s == 0.0) throw new java.lang.ArithmeticException("null multivector");
        else return this.gp(1 / Math.sqrt(Math.abs(s)));
    }*/
    /**
     * Unit under reverse norm.
     * 
     * nach Kleppe
     * normalize = {
     *   _P(1)/(sqrt(abs(_P(1)*_P(1)~)))
     * }
     *
     * @throws java.lang.ArithmeticException if multivector is null-vector
     * @param M metric
     * @return unit under 'reverse' norm (this / sqrt(abs(this.reverse(this))))
     */
    public Multivector unit_r(Metric M) {
        double s = scp(reverse(), M);
        if (s == 0.0) throw new java.lang.ArithmeticException("null multivector");
        else return this.gp(1 / Math.sqrt(Math.abs(s)));
    }

    protected double norm_e() {
        double s = scp(reverse());
        if (s < 0.0) return 0.0; // avoid FP round off causing negative 's'
        else return Math.sqrt(s);
    }
    
    public double norm_e(Metric m) {
        double s = scp(reverse(),m);
        if (s < 0.0) return 0.0; // avoid FP round off causing negative 's'
        else return Math.sqrt(s);
    }
    protected double norm_e2() {
        double s = scp(reverse());
        if (s < 0.0) return 0.0; // avoid FP round off causing negative 's'
        return s;
    }
    public double norm_e2(Metric m) {
        double s = scp(reverse(), m);
        if (s < 0.0) return 0.0; // avoid FP round off causing negative 's'
        return s;
    }
    
    /** 
     * @return true if this is really 0.0 
     */
    public boolean isNull() {
        simplify();
        return (blades.isEmpty());
    }

    /**
     * @param epsilon 
     *  
     * @return true if norm_e2 < epsilon * epsilon
     */
    public boolean isNull(double epsilon) {
        double s = norm_e2();
        return (s < epsilon * epsilon);
    }

    /** 
     * Is it a scalar.
     * 
     * @return true is this is a scalar (0.0 is also a scalar) 
     */
    public boolean isScalar() {
        if (isNull()) {
            return true;
        } else if (blades.size() == 1) {
            ScaledBasisBlade b = blades.get(0);
            return (b.bitmap == 0);
        } else {
            return false;
        }
    }

    /** 
     * Reverse of the multivector.
     * 
     * @return reverse of this 
     */
    public Multivector reverse() {
        List<ScaledBasisBlade> result = new ArrayList(blades.size());

        // loop over all basis lades, reverse them, add to result
        for (int i = 0; i < blades.size(); i++)
            result.add((blades.get(i)).reverse());

        return new Multivector(result);
    }

    /**
     * Grade inversion (is called involution in Dorst2007?) of the multivector.
     * 
     * This is also called "involution" in Dorst2007.
     * 
     * @return grade inversion/involution of this 
     */
    public Multivector gradeInversion() {
        List<ScaledBasisBlade> result = new ArrayList(blades.size());
        for (int i = 0; i < blades.size(); i++)
            result.add((blades.get(i)).gradeInversion());
        return new Multivector(result);
    }
    
    /**
     * Determines the grade of the basis blade with the maximum grade.
     * 
     * @return grade of the blade with the maximum grade
     */
    public int maxgrade(){
        int result = 0;
        for (int i = 0; i < blades.size(); i++){
            int grade = blades.get(i).grade();
            if (grade > result){
                result = grade;
            }
        }
        return result;
    }

    /**
     * @return clifford conjugate of this 
     * 
     * The conjugate operation is an involution, which means that twice application 
     * results in the original object.
     * 
     * The conjugate is equivalent to the reverse if all blades does not multiplay to -1.
     * This is not the case for CGA.
     */
    public Multivector cliffordConjugate() {
        List<ScaledBasisBlade> result = new ArrayList(blades.size());
        for (int i = 0; i < blades.size(); i++)
                result.add((blades.get(i)).cliffordConjugate());
        return new Multivector(result);
    }

    /**
     * Extracts grade 'g' from this multivector.
     * 
     * @param g
     * @return a new multivector of grade 'g'
     */
    public Multivector extractGrade(int g) {
        return extractGrade(new int[]{g});
    }
     
    /**
     * Extracts grade(s) 'G' from this multivector.
     * 
     * @param G
     * @return a new multivector of grade(s) 'G'
     */
    public Multivector extractGrade(int[] G) {
        return new Multivector(extractBlades(G));
    }

   /**
    * Extracts grade(s) 'G' from this multivector.
    * 
    * @param G
    * @return a new multivector of grade(s) 'G'
    */
    public List<ScaledBasisBlade> extractBlades(int[] G) {
        // what is the maximum grade to be extracted?
        int maxG = 0;
        for (int i = 0; i < G.length; i++)
            if (G[i] > maxG) maxG = G[i];

        // create boolean array of what grade to keep
        boolean[] keep = new boolean[maxG + 1]; // +1 to include the scalar
        for (int i = 0; i < G.length; i++)
            keep[G[i]] = true;

        // extract the grade, store in result:
        List<ScaledBasisBlade> result = new ArrayList();
        for (int i = 0; i < blades.size(); i++) {
            ScaledBasisBlade b = blades.get(i);
            int g = b.grade();
            if (g <= maxG && keep[g]) {
                result.add(b.copy());
            }
        }

        return result;
    }
    
    /**
     * Determination of the dual multivector.
     * 
     * @param dim
     * @return dual multivector
     */
    protected Multivector dual(int dim) {
        Multivector I = new Multivector(new ScaledBasisBlade((1 << dim)-1, 1.0));
        return ip(I.versorInverse(), LEFT_CONTRACTION);
    }

    /**
     * Determination of the dual multivector.
     * 
     * Keep in mind: Undualize is not the same as dual. It depends on the dimension
     * if there occures a sign.
     * 
     * Keep in mind: Depending on the dimension this can produce a different sign
     * than the dual()-method.
     * 
     * Dorst2007 page 80
     * @param M metric
     * @return dual multivector
     */
    public Multivector dual(Metric M) {
        Multivector I = new Multivector(new ScaledBasisBlade((1 << M.getEigenMetric().length)-1, 1.0));
        return ip(I.versorInverse(), M, LEFT_CONTRACTION);
    }
    
    /**
     * If the multivector is defined in dual representation the undual can be
     * determind.
     * 
     * Keep in mind: Depending on the dimension this can produce a different sign
     * than the dual()-method.
     * 
     * Dorst2007 page 80
     * 
     * FIXME
     * stimmt das überhaupt?
     * 
     * @param M metric
     * @return mulitvector in inner product null space representation if the original
     * multivector is defined in outer product null space representation.
     */
    public Multivector undual(Metric M){
         Multivector I = new Multivector(new ScaledBasisBlade((1 << M.getEigenMetric().length)-1, 1.0));
        return ip(I, M, LEFT_CONTRACTION);
    }

    /**
     * Determination of the dual multivector.
     * 
     * @param m metric
     * @return dual multivector
     */
    protected Multivector dual(double[] m) {
        Multivector I = new Multivector(new ScaledBasisBlade((1 << m.length)-1, 1.0));
        return ip(I.versorInverse(), m, LEFT_CONTRACTION);
    }

    /** 
     * Get the scalar part of the multivector. 
     * 
     * @return scalar part of 'this'
     */
    public double scalarPart() {
        double s = 0.0;
        for (int i = 0; i < blades.size(); i++) {
            ScaledBasisBlade b = blades.get(i);
            if (b.bitmap == 0) s += b.scale;
        }
        return s;
    }

    /**
     * Inverse of the multivector if it is a versor.
     * 
     * @return inverse of this (assuming it is a versor, no metric is used to 
     * check for invertability, no check to be a versor is made!)
     * @throws java.lang.ArithmeticException if versor is not invertible.
     */
    protected Multivector versorInverse() {
        Multivector R = reverse();
        double s = scp(R);
        if (s == 0.0) throw new java.lang.ArithmeticException("non-invertible multivector");
        return R.gp(1d / s);
    }

    /**
     * A versor is a multivector that can be expressed as the geometric product 
     * of a number of non-null 1-vectors. 
     * 
     * A sum of two versors does not in general result in a versor!<p>
     * 
     * @param M metric
     * @return inverse of this (assuming, it is a versor, no check is made!)
     * @throws java.lang.ArithmeticException if the multivector is not invertable
     */
    public Multivector versorInverse(Metric M) {
        Multivector R = reverse();
        double s = scp(R, M);
        if (s == 0.0) throw new java.lang.ArithmeticException("non-invertible multivector");
        return R.gp(1.0 / s);
    }

    /**
     * Can throw java.lang.ArithmeticException if versor is not invertible.
     * 
     * @param m
     * @return inverse of this (assuming it is a versor, no check is made!)
     * @throws java.lang.ArithmeticException if the multivector is not invertable
     */
    protected Multivector versorInverse(double[] m) {
        Multivector R = reverse();
        double s = scp(R, m);
        if (s == 0.0) throw new java.lang.ArithmeticException("non-invertible multivector");
        return R.gp(1.0 / s);
    }

    /**
     * Determine general inverse of this multivector.
     * 
     * @param metric
     * @return inverse of arbitrary multivector.
     * @throws java.lang.ArithmeticException if multivector is not invertible.
     */
    public Multivector generalInverse(Object metric) {
        int dim = spaceDim();

        cern.colt.matrix.DoubleMatrix2D M =
        cern.colt.matrix.DoubleFactory2D.dense.make(1 << dim, 1 << dim);

        // create all unit basis gblades for 'dim'
        ScaledBasisBlade[] B = new ScaledBasisBlade[1 << dim];
        for (int i = 0; i < (1 << dim); i++)
            B[i] = new ScaledBasisBlade(i);

        // construct a matrix 'M' such that matrix multiplication of 'M' with
        // the coordinates of another multivector 'x' (stored in a vector)
        // would result in the geometric product of 'M' and 'x'
        for (int i = 0; i < blades.size(); i++) {
            ScaledBasisBlade b = blades.get(i);
            for (int j = 0; j < (1 << dim); j++) {
                if (metric == null){
                    addToMatrix(M, b, B[j], ScaledBasisBlade.gp(b, B[j]));
                } else if (metric instanceof Metric metric1){
                    addToMatrix(M, b, B[j], ScaledBasisBlade.gp(b, B[j], metric1));
                } else if (metric instanceof double[] ds) {
                    addToMatrix(M, b, B[j], ScaledBasisBlade.gp(b, B[j], ds));
                }
            }
        }

        // try to invert matrix (can fail, then we throw an exception)
        cern.colt.matrix.DoubleMatrix2D IM = null;
        try {
            IM = cern.colt.matrix.linalg.Algebra.DEFAULT.inverse(M);
        } catch (Exception ex) {
            throw new java.lang.ArithmeticException("Multivector is not invertible");
        }

        // reconstruct multivector from first column of matrix
        List<ScaledBasisBlade> result = new ArrayList();
        for (int j = 0; j < (1 << dim); j++) {
            double v = IM.getQuick(j, 0);
            if (v != 0.0) {
                B[j].scale = v;
                result.add(B[j]);
            }
        }
        return new Multivector(result);
    }

    protected static void addToMatrix(cern.colt.matrix.DoubleMatrix2D M,
        ScaledBasisBlade alpha, ScaledBasisBlade beta, ScaledBasisBlade gamma) {
        // add gamma.scale to matrix entry M[gamma.bitmap, beta.bitmap]
        double v = M.getQuick(gamma.bitmap, beta.bitmap);
        M.setQuick(gamma.bitmap, beta.bitmap, v + gamma.scale);
    }

    protected static void addToMatrix(cern.colt.matrix.DoubleMatrix2D M,
                            ScaledBasisBlade alpha, ScaledBasisBlade beta, List<ScaledBasisBlade> gamma) {
        for (int i = 0; i < gamma.size(); i++)
            addToMatrix(M, alpha, beta, gamma.get(i));
    }

    /** 
     * Simplify the multivector representation.
     * 
     * @return simplification of this multivector (the same Multivector, 
     * but gblades array can be changed) 
     */
    public Multivector simplify() {
        simplify(blades);
        return this;
    }

    /** 
     * Get largest coordinate value of all gblades.
     * 
     * @return abs(largest coordinate) of all gblades.
     */
    public double largestCoordinate() {
        simplify();
        double lc = 0.0;
        for (int i = 0; i < blades.size(); i++) {
            ScaledBasisBlade b = blades.get(i);
            lc = Math.max(Math.abs(b.scale), lc);
        }
        return lc;
    }

    /** 
     * Get Basisblade with the largest coordinate.
     * 
     * @return abs(largest ScaledBasisBlade) of 'this' 
     */
    public ScaledBasisBlade largestBasisBlade() {
        simplify();
        ScaledBasisBlade bestBlade = null;
        double bestScale = -1.0;
        for (int i = 0; i < blades.size(); i++) {
            ScaledBasisBlade b = blades.get(i);
            if (Math.abs(b.scale) > bestScale) {
                bestScale = Math.abs(b.scale);
                bestBlade = b;
            }
        }
        return bestBlade;
    }

    /** 
     * @return the grade of this if homogeneous, -1 otherwise.
     * 0 is return for null Multivectors.
     */
    public int grade() {
        int g = -1;
        for (int i = 0; i < blades.size(); i++) {
            ScaledBasisBlade b = blades.get(i);
            if (g < 0) g = b.grade();
            else if (g != b.grade()) return -1;
        }
        return (g < 0) ? 0 : g;
    }


    /** 
     * @return bitmap of grades that are in use in 'this'
     */
    public int gradeUsage() {
        int gu = 0;
        for (int i = 0; i < blades.size(); i++) {
            ScaledBasisBlade b = blades.get(i);
            gu |= 1 << b.grade();
        }
        return gu;
    }

    /** 
     * @return index of highest grade in use in 'this'
     */
    public int topGradeIndex() {
        int maxG = 0;
        for (int i = 0; i < blades.size(); i++) {
            ScaledBasisBlade b = blades.get(i);
            maxG = Math.max(b.grade(), maxG);
        }
        return maxG;
    }

    /** 
     * @return the largest grade part of this 
     */
    public Multivector largestGradePart() {
        simplify();
        Multivector maxGP = null;
        double maxNorm = -1.0;
        int gu = gradeUsage();
        for (int i = 0; i <= topGradeIndex(); i++) {
            if ((gu & (1 << i)) == 0) continue;
            Multivector GP = extractGrade(i);
            double n = GP.norm_e();
            if (n > maxNorm) {
                maxGP = GP;
                maxNorm = n;
            }
        }
        return (maxGP == null) ? new Multivector() : maxGP;
    }

    /** 
     * Estimated space dimension (based on the used basisBlades).
     * 
     * FIXME
     * Eventuell im Konstruktor des Multivector spacedim als Konstante einführen.
     * 
     * @return dimension of space this blade (apparently) lives in 
     */
    protected int spaceDim() {
        int maxD = 0;
        for (int i = 0; i < blades.size(); i++) {
            ScaledBasisBlade b = blades.get(i);
            maxD = Math.max(Bits.highestOneBit(b.bitmap), maxD);
        }
        return maxD+1;
    }

    /**
     * Currently removes all basis gblades with |scale| less than epsilon.
     *
     * Old version did this:
     * Removes basis gblades with whose |scale| is less than <code>epsilon * maxMag</code> where
     * maxMag is the |scale| of the largest basis blade.
     *
     * @param epsilon
     * @return 'Compressed' version of this (the same Multivector, but gblades array can be changed)
     */
    public Multivector compress(double epsilon) {
        simplify();
        // find maximum magnitude:
        double maxMag = 0.0;
        for (int i = 0; i < blades.size(); i++) {
            ScaledBasisBlade b = blades.get(i);
            maxMag = Math.max(Math.abs(b.scale), maxMag);
        }
        if (maxMag == 0.0) {
            blades.clear();
        } else {
            // premultiply maxMag
            maxMag = epsilon; // used to read *=
            // remove basis gblades with too low scale
            for (Iterator<ScaledBasisBlade> I = blades.iterator(); I.hasNext(); ) {
                ScaledBasisBlade b = I.next();
                if (Math.abs(b.scale) < maxMag) I.remove();
            }
        }
        return this;
    }

    /** 
     * shortcut to compress(1e-13)
     * 
     * @return  
     */
    public Multivector compress() {
	return compress(1e-13);
    }
    
    public List<ScaledBasisBlade> getBlades() {
        return new ArrayList<>(blades);
    }
   
    /** sorts by bitmap only */
    private static class BladesComperator implements java.util.Comparator {
        @Override
        public int compare(Object o1, Object o2) {
            ScaledBasisBlade b1 = (ScaledBasisBlade)o1;
            ScaledBasisBlade b2 = (ScaledBasisBlade)o2;
            if (b1.bitmap < b2.bitmap) return -1;
            else if (b1.bitmap > b2.bitmap) return 1;
            else return 0;
        }
    }

    /**
     * Simplifies list of basis gblades; List is modified in the process 
     * 
     * @param L
     * @return simplified list of gblades
     */
    protected static List<ScaledBasisBlade> simplify(List<ScaledBasisBlade> L) {
        Collections.sort(L, new BladesComperator()); // sort by bitmap only
        ScaledBasisBlade prevBlade = null;
        boolean removeNullBlades = false;
        for (Iterator<ScaledBasisBlade> I = L.iterator(); I.hasNext(); ) {
            ScaledBasisBlade curBlade = I.next();
            if (curBlade.scale == 0.0) {
                I.remove();
                prevBlade = null;
            } else if ((prevBlade != null) && (prevBlade.bitmap == curBlade.bitmap)) {
                prevBlade.scale += curBlade.scale;
                I.remove();
            } else {
                if ((prevBlade != null) && (prevBlade.scale == 0.0))
                        removeNullBlades = true;
                prevBlade = curBlade;
            }
        }

        if (removeNullBlades) {
            for (Iterator<ScaledBasisBlade> I = L.iterator(); I.hasNext(); ) {
                ScaledBasisBlade curBlade = I.next();
                if (curBlade.scale == 0.0) I.remove();
            }
        }
        return L;
    }

    /** 
     * Replaces all BasisBlades found in array L
     * 
     * @param L
     * @return  
     */
    protected static List<ScaledBasisBlade> cloneArrayElements(List<ScaledBasisBlade> L) {
        for (int i = 0; i < L.size(); i++) {
            L.set(i, (ScaledBasisBlade) (L.get(i)).clone());
        }
        return L;
    }

    /** 
     * sorts the blade in 'gblades' based on bitmap only 
     */
    protected void sortBlades() {
        if (!bladesSorted) {
            Collections.sort(blades, new BladesComperator());
            bladesSorted = true;
        }
    }

    /** 
     * For internal use; M can be null, Metric or do.
     * 
     * @param x
     * @param M
     * @return  
     */
    protected Multivector gp(Multivector x, Object M) {
        if (M == null) {
            return gp(x);
        } else if (M instanceof Metric metric) {
            return gp(x, metric);
        } else {
            return gp(x, (double[]) M);
        }
    }

    /** 
     * For internal use only.
     * 
     * @param x
     * @param M null or double[]
     * @return  
     */
    protected double scp(Multivector x, Object M) {
        if (M == null) {
            return scp(x);
        } else if (M instanceof Metric metric){
            return scp(x, metric);
        } else {
            return scp(x, (double[])M);
        }
    }

    /** 
     * For internal use.
     * 
     * @param M null or Metric double[]
     * @return  
     */
    protected Multivector _versorInverse(Object M) {
        if (M == null) {
            return versorInverse();
        } else if (M instanceof Metric metric) {
            return versorInverse(metric);
        } else {
            return versorInverse((double[])M);
        }
    }
    
    /** list of basis gblades */
    protected List<ScaledBasisBlade> blades = null;

    /** when true, the gblades have been sorted on bitmap */
    protected boolean bladesSorted = false;
} 
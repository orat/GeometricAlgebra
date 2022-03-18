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

/*
 * MultivectorInfo.java
 *
 * Created on October 19, 2005, 12:14 PM
 *
 * Copyright 2005-2007 Daniel Fontijne, University of Amsterdam
 * fontijne@science.uva.nl
 *
 */
package de.orat.math.ga.basis;

import de.orat.math.ga.util.Bits;

/**
 * A class that computes whether a multivector is a blade or a versor,
 * and if so, what is the grade and parity of the blade/versor.
 *
 * <p>If you just want to know if something is a blade, don't provide
 * any metric, since non-invertible blades are blades too.
 *
 * @author  fontijne
 */
public final class MultivectorType {
    protected final static int _VERSOR = 2;
    protected final static int _BLADE = 4;

    public final static int MULTIVECTOR = 1;
    public final static int VERSOR = MULTIVECTOR | _VERSOR; // a versor is a multivector too
    public final static int BLADE = VERSOR | _BLADE; // a blade is a versor too


    /** true if multivector is zero */
    protected boolean zero;

    /** MULTIVECTOR, VERSOR, or BLADE */
    protected int type;
    /** top grade occupied by the multivector */
    protected int grade;
    /** grade usage (bitmap containing what grade parts are present)  */
    protected int gradeUsage;
    /** parity (1 = odd, 0 = even , -1 is none)*/
    protected int parity;

    /** 
     * Creates a new instance of MultivectorInfo.
     * 
     * @param A 
     */
    public MultivectorType(Multivector A) {
	init(A, null, 0.0);
    }

    /** 
     * Creates a new instance of MultivectorInfo.
     * 
     * @param A
     * @param epsilon 
     */
    public MultivectorType(Multivector A, double epsilon) {
	init(A, null, epsilon);
    }

    /** 
     * Creates a new instance of MultivectorType.
     * 
     * @param A
     * @param metric */
    public MultivectorType(Multivector A, Object metric) {
	init(A, metric, 0.0);
    }

    /** 
     * Creates a new instance of MultivectorType.
     * 
     * @param A
     * @param metric
     * @param epsilon */
    public MultivectorType(Multivector A, Object metric, double epsilon) {
	init(A, metric, epsilon);
    }

    /**
     * Called by constructor to initialize the instance.
     * 
     * Performs the actual test if 'A' is a versor or blade under 'metric'
     * 
     * @param A
     * @param metric
     * @param epsilon
     */
    protected void init(Multivector A, Object metric, double epsilon) {
	// first of all determine grade usage & parity
	Multivector cA = ((epsilon == 0.0) ? A : A.compress(epsilon));
	grade = cA.topGradeIndex();
	gradeUsage = cA.gradeUsage();
	{
	    final int[] cnt = new int[]{0, 0}; // nb even, odd grade parts in use
	    int cntIdx = 0;
	    int gu = gradeUsage;
	    while (gu != 0) {
		if ((gu & 1) != 0)
		    cnt[cntIdx]++;
		cntIdx ^= 1;
		gu >>>= 1;
	    }

	    if ((cnt[0] == 0) && (cnt[1] == 0)) {
		// multivector = zero blade, case closed
		zero = true;
		type = BLADE;
		parity = 0;
		grade = -1;
	        return;
	    } else {
		zero = false;
		if ((cnt[0] != 0) && (cnt[1] != 0)) {
		    // multivector = multivector, case closed
		    type = MULTIVECTOR;
		    parity = -1;
		    return;
		} else // more work to do, but parity is known:
		    parity = (cnt[1] != 0) ? 1 : 0;
	    }
	}

	// determine if versor or blade:

	// get versor inverse, grade inversion:
	Multivector Avi;
	try {
	    Avi = A._versorInverse(metric);
	} catch (java.lang.ArithmeticException ex) {
	    type = MULTIVECTOR;
	    return;
	}
	Multivector Agi = A.gradeInversion();

	// check if Agi Avi is a scalar:
	if (Agi._gp(Avi, metric).compress(epsilon).gradeUsage() != 1) {
	    // multivector = multivector, case closed
	    type = MULTIVECTOR;
	    return;
	}
	// check if Agi Avi == Avi Agi
	if (!Agi._gp(Avi, metric).subtract(Avi._gp(Agi, metric)).compress(epsilon).isNull()) {
	    // multivector = multivector, case closed
	    type = MULTIVECTOR;
	    return;
	}

	// check if grade preserving for all basis vectors:
	int dim = A.spaceDim();
	for (int i = 0; i < dim; i++) {
	    Multivector ei = new Multivector(new ScaledBasisBlade(1 << i));
	    if (Avi._gp(ei, metric)._gp(Agi, metric).compress().gradeUsage() != 2) {
		// multivector = multivector, case closed
		type = MULTIVECTOR;
		return;
	    }

	}

	// if homogeneous: blade
	if (Bits.bitCount(gradeUsage) == 1)
	    type = BLADE;
	else type = VERSOR;
    }

    /** 
     * @return MULTIVECTOR, VERSOR, or BLADE 
     */
    public int getType() {
	return type;
    }

    /** 
     * @return true if multivector is zero 
     */
    public boolean isZero() {
	return zero;
    }

    /** 
     * @return true if multivector is a versor 
     */
    public boolean isVersor() {
	return ((getType() & _VERSOR) != 0);
    }

    /** 
     * Is blade.
     * 
     * @return true if multivector is a blade 
     */
    public boolean isBlade() {
	return ((getType() & _BLADE) != 0);
    }

    /** 
     * Get top grade index.
     * 
     * @return grade of the multivector 
     */
    public int getTopGradeIndex() {
	return grade;
    }

    /** 
     * Get grade usage.
     * 
     * @return grade usage of the multivector 
     */
    public int getGradeUsage() {
	return gradeUsage;
    }

    /** 
     * @return true if multivector is even versor 
     */
    public boolean isEven() {
	return parity == 0;
    }
    public boolean isOdd() {
	return parity == 1;
    }

    /** 
     * parity 
     * 
     * @return 1 = odd, 0 = even , -1 is none
     */
    public int getParity() {
	return parity;
    }

    @Override
    public String toString() {
	return "Multivector type[" +
	 "type = " + (isBlade() ? "blade" : (isVersor() ? "versor" : "multivector")) +
	 ", parity = " + (isEven() ? "even" : (isOdd() ? "odd" : "none")) +
	 ", top grade index = " + getTopGradeIndex() +
	 ", grade usage = " + getGradeUsage() + "]";
    }
}

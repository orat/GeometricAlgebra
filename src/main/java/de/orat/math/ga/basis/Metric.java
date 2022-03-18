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
 * Metric.java
 *
 * Created on February 1, 2005, 12:20 PM
 *
 * Copyright 2005-2007 Daniel Fontijne, University of Amsterdam
 * fontijne@science.uva.nl
 *
 */

package de.orat.math.ga.basis;

import java.util.ArrayList;
import cern.colt.matrix.*;
import cern.colt.matrix.linalg.*;
import de.orat.math.ga.metric.MetricException;
import java.util.List;

/**
 * A class for representing totally arbitrary metric. 
 * 
 * Metric is given
 * in the form of a (symmetric) matrix. The eigenvectors of the matrix
 * are computed since these are used for computing (geometric) products
 * in this metric.
 *
 * @author  fontijne
 */
public class Metric implements de.orat.math.ga.metric.Metric {
    public static void main(String args[]) {
        double[][] m = new double[][]{
            {1.0, 0.0, 0.0},
            {0.0 ,0.0, -1.0},
            {0.0 , -1.0, 0.0}
        };

        try {
            Metric M = new Metric(m);
        } catch (MetricException E) {
        }

    }

    /** Creates a new instance of Metric
     * @param m
     * @throws de.orat.math.ga.impl.metric.MetricException */
    public Metric(double[][] m) throws MetricException {
        this(DoubleFactory2D.dense.make(m));
    }

    /** 
     * Creates a new instance of Metric.
     * 
     * @param m
     * @throws de.orat.math.ga.impl.metric.MetricException 
     */
    public Metric(DoubleMatrix2D m) throws MetricException {
        m_matrix = m.copy();
        if (!Property.ZERO.isSymmetric(m_matrix))
            throw new MetricException("The metric matrix must be symmetric");

        // compute eigen value decomposition
        m_eig = new EigenvalueDecomposition(m_matrix);

        m_invEigMatrix = Algebra.ZERO.transpose(m_eig.getV());

        m_eigenMetric = new double[m_eig.getRealEigenvalues().size()];
        for (int i = 0; i < m_eigenMetric.length; i++)
            m_eigenMetric[i] = m_eig.getRealEigenvalues().get(i);


		m_isDiagonal = cern.colt.matrix.linalg.Property.ZERO.isDiagonal(m_matrix);
		if (!m_isDiagonal) {
			m_isEuclidean = m_isAntiEuclidean = false;
		}
		else {
			m_isEuclidean = m_isAntiEuclidean = true;
			for (int i = 0; i < m_matrix.columns(); i++) {
			if (m_matrix.getQuick(i, i) != 1.0)
				m_isEuclidean = false;
			if (m_matrix.getQuick(i, i) != -1.0)
				m_isAntiEuclidean = false;
			}
		}
    }

    /** transforms a to the eigen basis (a must be on metric basis) */
    public List<ScaledBasisBlade> toEigenbasis(ScaledBasisBlade a) {
        return transform(a, m_invEigMatrix);
    }

    /** 
     * transforms A to the eigen basis (A must be on metric basis)
     * 
     * @param a
     * @return  
     */
    public List<ScaledBasisBlade> toEigenbasis(List<ScaledBasisBlade> a) {
        List<ScaledBasisBlade> result = new ArrayList();
        for (int i = 0; i < a.size(); i++) {
            List tmp = toEigenbasis(a.get(i));
            result.addAll(tmp);
        }
        return ScaledBasisBlade.simplify(result);
    }

    /** 
     * transforms a to the metric basis (A must be on eigenbasis).
     * 
     * @param a
     * @return  
     */
    public List<ScaledBasisBlade> toMetricBasis(ScaledBasisBlade a) {
        return transform(a, m_eig.getV());
    }

    /** transforms A to the metric basis (a must be on eigenbasis).
     * 
     * @param a
     * @return  */
    public List<ScaledBasisBlade> toMetricBasis(List<ScaledBasisBlade> a) {
        List<ScaledBasisBlade> result = new ArrayList();
        for (int i = 0; i < a.size(); i++) {
            List<ScaledBasisBlade> tmp = toMetricBasis(a.get(i));
            result.addAll(tmp);
        }
        return ScaledBasisBlade.simplify(result);
    }

    /** 
     * transforms a ScaledBasisBlade to a new basis.
     * 
     * @param a
     * @param M
     * @return  */
    protected List<ScaledBasisBlade> transform(ScaledBasisBlade a, DoubleMatrix2D M) {
        ArrayList A = new ArrayList();
        A.add(new ScaledBasisBlade(a.scale)); // start with just scalar;

        // for each 1 bit: convert to list of blades
        int i = 0;
        int b = a.bitmap;
        while (b != 0) {
            if ((b & 1) != 0) {
                // take column 'i' out of the matrix, wedge it to 'A'
                ArrayList tmp = new ArrayList();
                for (int j = 0; j < M.rows(); j++) {
                    if (M.get(j, i) != 0) {
                        double m = M.get(j, i);
                        for (int k = 0; k < A.size(); k++)
                            tmp.add(ScaledBasisBlade.op((ScaledBasisBlade)A.get(k), new ScaledBasisBlade((1 << j), m)));
                    }
                }
                A = tmp;
            }

            b >>= 1;
            i++;
        }
        return A;
    }

    /** 
     * returns metric when in 'eigenspace'. Do not modify the array! 
     */
    public double[] getEigenMetric() {
        return m_eigenMetric;
    }

    public double getBasisVectorIp(int idx) {
		return getBasisVectorIp(idx, idx);
    }

    public double getBasisVectorIp(int idx1, int idx2) {
		return m_matrix.get(idx1, idx2);
    }

    public cern.colt.matrix.DoubleMatrix2D getInnerProductMatrix() {
		return m_matrix;
    }

    public cern.colt.matrix.linalg.EigenvalueDecomposition getInnerProductMatrixEigenvalueDecomposition() {
		return m_eig;
    }

    public boolean isAntiEuclidean() {
		return m_isAntiEuclidean;
    }

    public boolean isEuclidean() {
		return m_isEuclidean;
    }

    public boolean isDiagonal() {
	return m_isDiagonal;
    }

    /** the metric (symmetric matrix) */
    private DoubleMatrix2D m_matrix;

    /** the eigenvectors matrix & eigenvalues of m_matrix */
    protected EigenvalueDecomposition m_eig;

    /** inverse of the eigenmatrix */
    protected DoubleMatrix2D m_invEigMatrix;

    protected double[] m_eigenMetric;

    protected boolean m_isDiagonal;
    protected boolean m_isEuclidean;
    protected boolean m_isAntiEuclidean;
}

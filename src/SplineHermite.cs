using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Landis.Extension.GiplDamm
{
    public static class SplineHermite
    {
        public static void SplineHermiteSet(int ndata, double[] tdata, double[] ydata, double[] ypdata, double[][] c)
        {
            //% SPLINE_HERMITE_SET sets up a piecewise cubic Hermite interpolant.
            //
            // Discussion:
            //
            // Once the array C is computed, then in the interval
            // (TDATA(I), TDATA(I + 1)), the interpolating Hermite polynomial
            //    is given by
            //
            // SVAL(TVAL) = C(1, I)
            // +(TVAL - TDATA(I)) * (C(2, I)
            // +(TVAL - TDATA(I)) * (C(3, I)
            // +(TVAL - TDATA(I)) * C(4, I) ) )
            //
            // Licensing:
            //
            // This code is distributed under the GNU LGPL license.
            //
            // Modified:
            //
            // 11 February 2004
            //
            // Author:
            //
            // John Burkardt
            //
            // Reference:
            //
            // Samuel Conte and Carl de Boor,
            // Algorithm CALCCF,
            // Elementary Numerical Analysis,
            // 1973, page 235.
            //
            // Parameters:
            //
            // Input, integer NDATA, the number of data points.
            // NDATA must be at least 2.
            //
            // Input, real TDATA(NDATA), the abscissas of the data points.
            // The entries of TDATA are assumed to be strictly increasing.
            //
            // Input, real Y(NDATA), YP(NDATA), the value of the
            // function and its derivative at TDATA(1:NDATA).
            //
            // Output, real C(4, NDATA), the coefficients of the Hermite polynomial.
            // C(1, 1:NDATA) = Y(1:NDATA) and C(2, 1:NDATA) = YP(1:NDATA).
            // C(3, 1:NDATA - 1) and C(4, 1:NDATA - 1) are the quadratic and cubic
            // coefficients.
            //

            Array.Copy(ydata, c[1], ndata + 1);
            Array.Copy(ypdata, c[2], ndata + 1);

            for (var i = 1; i < ndata; ++i)
            {
                var dt = tdata[i + 1] - tdata[i];
                var divdif1 = (c[1][i + 1] - c[1][i]) / dt;
                var divdif3 = c[2][i] + c[2][i + 1] - 2.0 * divdif1;
                c[3][i] = (divdif1 - c[2][i] - divdif3) / dt;
                c[4][i] = divdif3 / (dt * dt);
            }

            c[3][ndata] = 0.0;
            c[4][ndata] = 0.0;
        }

        public static void SplineHermitVal(int ndata, double[] tdata, double[][] c, double tval, out double sval, out double spval, out double sp2val)
        {
            //% SPLINE_HERMITE_VAL evaluates a piecewise cubic Hermite interpolant.
            //
            // Discussion:
            //
            // SPLINE_HERMITE_SET must be called first, to set up the
            // spline data from the raw function and derivative data.
            //
            // In the interval (TDATA(I), TDATA(I + 1)), the interpolating
            // Hermite polynomial is given by
            //
            // SVAL(TVAL) = C(1, I)
            // +(TVAL - TDATA(I)) * (C(2, I)
            // +(TVAL - TDATA(I)) * (C(3, I)
            // +(TVAL - TDATA(I)) * C(4, I) ) )
            //
            // and
            //
            // SVAL'(TVAL) =                    C(2,I)
            // +(TVAL - TDATA(I)) * (2 * C(3, I)
            // +(TVAL - TDATA(I)) * 3 * C(4, I) )
            //
            // Licensing:
            //
            // This code is distributed under the GNU LGPL license.
            //
            // Modified:
            //
            // 11 February 2004
            //
            // Author:
            //
            // John Burkardt
            //
            // Reference:
            //
            // Samuel Conte and Carl de Boor,
            // Algorithm PCUBIC,
            // Elementary Numerical Analysis,
            // 1973, page 234.
            //
            // Parameters:
            //
            // Input, integer NDATA, the number of data points.
            // NDATA is assumed to be at least 2.
            //
            // Input, real TDATA(NDATA), the abscissas of the data points.
            // The entries of TDATA are assumed to be strictly increasing.
            //
            // Input, real C(4, NDATA), the coefficient data computed by
            // SPLINE_HERMITE_SET.
            //
            // Input, real TVAL, the point where the interpolant is to
            // be evaluated.
            //
            // Output, real SVAL, SPVAL, the value of the interpolant
            // and its derivative at TVAL.

            var left = 0;
            var right = 0;

            // Find the interval[TDATA(LEFT), TDATA(RIGHT)] that contains
            // or is nearest to TVAL.

            Bracket(ndata, tdata, tval, ref left, ref right);

            // Evaluate the cubic polynomial.

            var dt = tval - tdata[left];
            sval = c[1][left] + dt * (c[2][left] + dt * (c[3][left] + dt * c[4][left]));
            spval = c[2][left] + dt * (2.0 * c[3][left] + dt * 3.0 * c[4][left]);
            sp2val = 2.0 * c[3][left] + dt * 6.0 * c[4][left];
        }

        private static void Bracket(int n, double[] x, double xval , ref int left, ref int right)
        {
            // %% R8VEC_BRACKET searches a sorted array for successive brackets of a value.
            // %
            // %  Discussion:
            // %
            // %    If the values in the vector are thought of as defining intervals
            // %    on the real line, then this routine searches for the interval
            // %    nearest to or containing the given value.
            // %
            // %  Licensing:
            // %
            // %    This code is distributed under the GNU LGPL license.
            // %
            // %  Modified:
            // %
            // %    24 January 2004
            // %
            // %  Author:
            // %
            // %    John Burkardt
            // %
            // %  Parameters:
            // %
            // %    Input, integer N, length of input array.
            // %
            // %    Input, real X(N), an array that has been sorted into ascending order.
            // %
            // %    Input, real XVAL, a value to be bracketed.
            // %
            // %    Output, integer LEFT, RIGHT, the results of the search.
            // %    Either:
            // %      XVAL < X(1), when LEFT = 1, RIGHT = 2;
            // %      XVAL > X(N), when LEFT = N-1, RIGHT = N;
            // %    or
            // %      X(LEFT) <= XVAL <= X(RIGHT).
            // %

            for (var i = 2; i <= n - 1; ++i)
            {
                if (xval < x[i])
                {
                    left = i - 1;
                    right = i;
                    goto label1;
                }
            }
            left = n - 1;
            right = n;
            label1:;
        }
    }
}

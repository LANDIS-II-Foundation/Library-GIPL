using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Landis.Extension.GiplDamm
{
    public partial class GiplDamm
    {
        private void Stemperature(IndexedArray<double> soil_temperature, IndexedArray<double> liquid_fraction, ref double time, double time_end, ref double time_step, ref int gc, double[] airTemp)
        {
            Stemperature_atomic(soil_temperature, liquid_fraction, ref time_step, ref time, time_end, ref gc, airTemp);
            var old_time_step = time_step;
            var dt_temp = time_end - time;

            while (1.0e-10 < dt_temp)
            {
                Stemperature_atomic(soil_temperature, liquid_fraction, ref dt_temp, ref time, time_end, ref gc, airTemp);
                dt_temp = time_end - time;
            }
            time_step = old_time_step;
        }

        private void Stemperature_atomic(IndexedArray<double> soil_temperature, IndexedArray<double> liquid_fraction, ref double dt, ref double time, double time_end, ref int gc, double[] airTemp)
        {
            // INCLUDE 'pwbm.soil_temperature.h'
            var it = 0;
            var flag = 0;
            var t = new IndexedArray<double>(_numericsModule.Sx, _numericsModule.Lx + 1);
            var wc = 0.0;
            var wp = 0.0;
            var dtemperature = 0.0;
            var capp_water = 0.0;
            var capp = 0.0;
            var delta = 0.0;
            var c1 = 0.0;
            var b1 = 0.0;
            var d1 = 0.0;
            var r1 = 0.0;
            var l1 = 0.0;
            var tl = new IndexedArray<double>(_numericsModule.Sx, _numericsModule.Lx + 1);
            var tr = new IndexedArray<double>(_numericsModule.Sx, _numericsModule.Lx + 1);
            var atr = new IndexedArray<double>(_numericsModule.Sx, _numericsModule.Lx + 1);
            var btr = new IndexedArray<double>(_numericsModule.Sx, _numericsModule.Lx + 1);
            var ctr = new IndexedArray<double>(_numericsModule.Sx, _numericsModule.Lx + 1);
            var atrnr = new IndexedArray<double>(_numericsModule.Sx, _numericsModule.Lx + 1);
            var btrnr = new IndexedArray<double>(_numericsModule.Sx, _numericsModule.Lx + 1);
            var ctrnr = new IndexedArray<double>(_numericsModule.Sx, _numericsModule.Lx + 1);
            var dtr = new IndexedArray<double>(_numericsModule.Sx, _numericsModule.Lx + 1);
            var dtrnr = new IndexedArray<double>(_numericsModule.Sx, _numericsModule.Lx + 1);
            // Heat capacity
            var c = new IndexedArray<double>(_numericsModule.Sx, _numericsModule.Lx + 1);
            var dc_dt = new IndexedArray<double>(_numericsModule.Sx, _numericsModule.Lx + 1);
            // Enthalpy due to water
            var h = new IndexedArray<double>(_numericsModule.Sx, _numericsModule.Lx + 1);
            var dh_dt = new IndexedArray<double>(_numericsModule.Sx, _numericsModule.Lx + 1);
            var d2h_d2t = new IndexedArray<double>(_numericsModule.Sx, _numericsModule.Lx + 1);
            // Thermal conductivity
            var k = new IndexedArray<double>(_numericsModule.Sx, _numericsModule.Lx + 1);
            var dk_dt = new IndexedArray<double>(_numericsModule.Sx, _numericsModule.Lx + 1);
            // Liquid water fraction
            var y = new IndexedArray<double>(_numericsModule.Sx, _numericsModule.Lx + 1);
            var dy_dt = new IndexedArray<double>(_numericsModule.Sx, _numericsModule.Lx + 1);
            var d2y_d2t = new IndexedArray<double>(_numericsModule.Sx, _numericsModule.Lx + 1);
            var surface_temperture = 0.0;

            // local transfer variables for indexed values within IndexedArrays that can't be passed as 'refs' or 'outs'
            double yi, dy_dti, d2y_d2ti, hi, dh_dti, d2h_d2ti, ci, dc_dti, ki, dk_dti;
            double liquid_fractioni;

            while (time + dt <= time_end)
            {
                for (var i = _numericsModule.Sn; i <= _numericsModule.Lx + 1; ++i)
                {
                    t[i] = soil_temperature[i];
                }
                flag = 1;
                it = 0;
                while (flag == 1)
                {
                    for (var i = _numericsModule.Sn; i <= 0; ++i)
                    {
                        Snow_th_properties(t[i], out yi, out dy_dti, out d2y_d2ti, out hi, out dh_dti, out d2h_d2ti, out ci, out dc_dti, out ki, out dk_dti, _numericsModule.LayerW[i], _numericsModule.LayerCf[i], _numericsModule.LayerLf[i]);
                        y[i] = yi; dy_dt[i] = dy_dti; d2y_d2t[i] = d2y_d2ti; h[i] = hi; dh_dt[i] = dh_dti; d2h_d2t[i] = d2h_d2ti; c[i] = ci; dc_dt[i] = dc_dti; k[i] = ki; dk_dt[i] = dk_dti;
                    }
                    for (var i = 2; i <= _numericsModule.Lx + 1; ++i)
                    {
                        Soil_th_properties(t[i], out yi, out dy_dti, out d2y_d2ti, out hi, out dh_dti, out d2h_d2ti, out ci, out dc_dti, out ki, out dk_dti, _numericsModule.LayerW[i], _numericsModule.LayerCf[i], _numericsModule.LayerLf[i], _numericsModule.LayerType[i]);
                        y[i] = yi; dy_dt[i] = dy_dti; d2y_d2t[i] = d2y_d2ti; h[i] = hi; dh_dt[i] = dh_dti; d2h_d2t[i] = d2h_d2ti; c[i] = ci; dc_dt[i] = dc_dti; k[i] = ki; dk_dt[i] = dk_dti;
                    }
                    // Snow model
                    for (var i = _numericsModule.Sn; i <= -1; ++i)
                    {
                        tl[i] = (k[i] + k[i + 1]) * dt * _numericsModule.Dx1[i] * 0.5;
                        btr[i] = tl[i];
                        ctr[i] = -tl[i];
                        atr[i + 1] = -tl[i];
                        btrnr[i] = tl[i];
                        ctrnr[i] = -tl[i];
                        atrnr[i + 1] = -tl[i];
                    }
                    btr[0] = 0.0;
                    btrnr[0] = 0.0;
                    // Flux
                    ctr[0] = 1.0;
                    ctrnr[0] = 1.0;
                    // Mortar
                    btr[1] = 0.0;
                    btrnr[1] = 0.0;
                    atr[1] = 1.0;
                    atrnr[1] = 1.0;
                    ctr[1] = -1.0;
                    ctrnr[1] = -1.0;
                    dtr[1] = 0.0;
                    // Flux
                    atr[2] = -1.0; atrnr[2] = -1.0;
                    atr[_numericsModule.Sn] = 0.0;
                    atrnr[_numericsModule.Sn] = 0.0;
                    // Soil
                    for (var i = 2; i <= _numericsModule.Lx; ++i)
                    {
                        tl[i] = (k[i] + k[i + 1]) * dt * _numericsModule.Dx1[i - 1] * 0.5;
                        tr[i] = (t[i] - t[i + 1]) * dt * _numericsModule.Dx1[i - 1] * 0.5;
                        btr[i] = tl[i];
                        ctr[i] = -tl[i];
                        atr[i + 1] = -tl[i];
                        btrnr[i] = tl[i] + tr[i] * dk_dt[i];
                        ctrnr[i] = -tl[i] + tr[i] * dk_dt[i + 1];
                        atrnr[i + 1] = -tl[i] - tr[i] * dk_dt[i];
                    }
                    btr[_numericsModule.Lx + 1] = 0.0;
                    btrnr[_numericsModule.Lx + 1] = 0.0;
                    ctr[_numericsModule.Lx + 1] = 0.0;
                    ctrnr[_numericsModule.Lx + 1] = 0.0;
                    // Snow
                    for (var i = _numericsModule.Sn; i <= -1; ++i)
                    {
                        btr[i + 1] = btr[i + 1] + tl[i];
                        btrnr[i + 1] = btrnr[i + 1] + tl[i];
                    }
                    // Soil
                    for (var i = 2; i <= _numericsModule.Lx; ++i)
                    {
                        btr[i + 1] = btr[i + 1] + tl[i];
                        btrnr[i + 1] = btrnr[i + 1] + tl[i] - tr[i] * dk_dt[i + 1];
                    }
                    // Snow
                    for (var i = _numericsModule.Sn; i <= 0; ++i)
                    {
                        btr[i] = btr[i] + c[i] * _numericsModule.Diag[i];
                        btrnr[i] = btrnr[i] + c[i] * _numericsModule.Diag[i];
                        dtr[i] = c[i] * _numericsModule.Diag[i] * soil_temperature[i];
                    }
                    // Soil
                    for (var i = 2; i <= _numericsModule.Lx + 1; ++i)
                    {
                        dtemperature = t[i] - soil_temperature[i];
                        // Second term in the apparent heat capacity due to freezing water
                        if (Math.Abs(dtemperature) < NumericsModule.Dlt)
                        {
                            // if temperature in the soil layer is the same we approximate it analytically
                            capp_water = dh_dt[i];
                            btrnr[i] = btrnr[i] + d2h_d2t[i] * _numericsModule.Diag[i - 1] * dtemperature;
                        }
                        else
                        {
                            // Enthalpy due to unfrozen water at the previous step
                            wp = _numericsModule.LayerW[i] * liquid_fraction[i] * NumericsModule.La * _numericsModule.Tscale;
                            // Enthalpy due to unfrozen water at the current step
                            wc = h[i];
                            // Otherwise we use an approximation
                            capp_water = (wc - wp) / dtemperature;
                            btrnr[i] = btrnr[i] - (capp_water - dh_dt[i]) * _numericsModule.Diag[i - 1];
                        }
                        // Apparent heat capacity at the current step
                        capp = c[i] + capp_water;
                        btr[i] = btr[i] + capp * _numericsModule.Diag[i - 1];
                        btrnr[i] = btrnr[i] + capp * _numericsModule.Diag[i - 1] + dc_dt[i] * _numericsModule.Diag[i - 1] * (t[i] - soil_temperature[i]);
                        dtr[i] = capp * _numericsModule.Diag[i - 1] * soil_temperature[i];
                    }
                    // DTR(SN)=DTR(SN)+flux*DT
                    r1 = SurfaceTemperature(time + dt, airTemp);
                    if (it < 1)
                    {
                        ctr[_numericsModule.Sn] = 0.0;
                        btr[_numericsModule.Sn] = 1.0;
                        dtr[_numericsModule.Sn] = r1;
                        Tridiag(_numericsModule.Sn, _numericsModule.Sx, _numericsModule.Lx + 1, atr, btr, ctr, dtr);
                        delta = 1.0;
                        delta = 0.0;
                        for (var i = _numericsModule.Sn; i <= _numericsModule.Lx + 1; ++i)
                        {
                            delta = Math.Max(Math.Abs(dtr[i] - t[i]), delta);
                        }
                        for (var i = _numericsModule.Sn; i <= _numericsModule.Lx + 1; ++i)
                        {
                            t[i] = dtr[i];
                        }
                    }
                    else
                    {
                        dtrnr[_numericsModule.Sn] = btr[_numericsModule.Sn] * t[_numericsModule.Sn] + ctr[_numericsModule.Sn] * t[_numericsModule.Sn + 1] - dtr[_numericsModule.Sn];
                        for (var i = _numericsModule.Sn + 1; i <= _numericsModule.Lx; ++i)
                        {
                            dtrnr[i] = atr[i] * t[i - 1] + btr[i] * t[i] + ctr[i] * t[i + 1] - dtr[i];
                        }
                        dtrnr[_numericsModule.Lx + 1] = atr[_numericsModule.Lx + 1] * t[_numericsModule.Lx] + btr[_numericsModule.Lx + 1] * t[_numericsModule.Lx + 1] - dtr[_numericsModule.Lx + 1];
                        ctrnr[_numericsModule.Sn] = 0.0;
                        btrnr[_numericsModule.Sn] = 1.0;
                        dtrnr[_numericsModule.Sn] = t[_numericsModule.Sn] - r1;
                        Tridiag(_numericsModule.Sn, _numericsModule.Sx, _numericsModule.Lx + 1, atrnr, btrnr, ctrnr, dtrnr);
                        delta = 0.0;
                        for (var i = _numericsModule.Sn; i <= _numericsModule.Lx + 1; ++i)
                        {
                            t[i] = t[i] - dtrnr[i];
                            delta = Math.Max(Math.Abs(dtrnr[i]), delta);
                        }
                    }
                    // WRITE(*,'(F10.6,F16.13 ,I3,F20.16,I4)') TIME+DT,DELTA,IT,DT,GC
                    if (delta < NumericsModule.MaxDelta)
                    {
                        time = time + dt;
                        gc = gc + 1;
                        if (gc == 5)
                        {
                            dt = dt * 2.0;
                            gc = 0;
                        }
                        dt = Math.Min(dt, NumericsModule.MaxDt);
                        it = 0;
                        flag = 0;
                    }
                    else
                    {
                        it = it + 1;
                    }
                    if (it == NumericsModule.NumMaxIter)
                    {
                        dt = dt / 2.0;
                        for (var i = _numericsModule.Sn; i <= _numericsModule.Lx + 1; ++i)
                        {
                            t[i] = soil_temperature[i];
                        }
                        gc = 0;
                        it = 0;
                    }
                }
                for (var i = _numericsModule.Sn; i <= 0; ++i)
                {
                    soil_temperature[i] = t[i];
                    Snow_th_properties(soil_temperature[i], out liquid_fractioni, out dy_dti, out d2y_d2ti, out hi, out dh_dti, out d2h_d2ti, out ci, out dc_dti, out ki, out dk_dti, _numericsModule.LayerW[i], _numericsModule.LayerCf[i], _numericsModule.LayerLf[i]);
                    liquid_fraction[i] = liquid_fractioni; dy_dt[i] = dy_dti; d2y_d2t[i] = d2y_d2ti; h[i] = hi; dh_dt[i] = dh_dti; d2h_d2t[i] = d2h_d2ti; c[i] = ci; dc_dt[i] = dc_dti; k[i] = ki; dk_dt[i] = dk_dti;
                }
                soil_temperature[1] = t[1];
                for (var i = 2; i <= _numericsModule.Lx + 1; ++i)
                {
                    soil_temperature[i] = t[i];
                    Soil_th_properties(soil_temperature[i], out liquid_fractioni, out dy_dti, out d2y_d2ti, out hi, out dh_dti, out d2h_d2ti, out ci, out dc_dti, out ki, out dk_dti, _numericsModule.LayerW[i], _numericsModule.LayerCf[i], _numericsModule.LayerLf[i], _numericsModule.LayerType[i]);
                    liquid_fraction[i] = liquid_fractioni; dy_dt[i] = dy_dti; d2y_d2t[i] = d2y_d2ti; h[i] = hi; dh_dt[i] = dh_dti; d2h_d2t[i] = d2h_d2ti; c[i] = ci; dc_dt[i] = dc_dti; k[i] = ki; dk_dt[i] = dk_dti;
                }
                for (var i = _numericsModule.Sx; i <= _numericsModule.Sn - 1; ++i)
                {
                    soil_temperature[i] = soil_temperature[_numericsModule.Sn];
                    liquid_fraction[i] = liquid_fraction[_numericsModule.Sn];
                }
            }
        }

        private void Snow_th_properties(double t, out double y, out double dy_dt, out double d2y_d2t, out double h, out double dh_dt, out double d2h_d2t, out double c, out double dc_dt, out double k, out double dk_dt, double w, double cf, double lf)
        {
            c = cf * _numericsModule.Tscale;
            k = lf;
            dk_dt = 0.0;
            dc_dt = 0.0;
            y = 1.0;
            dy_dt = 0.0;
            d2y_d2t = 0.0;
            h = 0.0;
            dh_dt = 0.0;
            d2h_d2t = 0.0;
        }

        private void Soil_th_properties(double t, out double y, out double dy_dt, out double d2y_d2t, out double h, out double dh_dt, out double d2h_d2t, out double c, out double dc_dt, out double k, out double dk_dt, double w, double cf, double lf, int shg)
        {
            Soilsat(t, out y, out dy_dt, out d2y_d2t, shg);
            c = (cf + w * y * (NumericsModule.Cw - NumericsModule.Ci)) * _numericsModule.Tscale;
            dc_dt = w * dy_dt * (NumericsModule.Cw - NumericsModule.Ci) * _numericsModule.Tscale;
            var kwki = Math.Exp(y * w * NumericsModule.Lkwki);
            k = lf * kwki;
            dk_dt = lf * kwki * w * NumericsModule.Lkwki * dy_dt;
            h = w * y * NumericsModule.La * _numericsModule.Tscale;
            dh_dt = w * dy_dt * NumericsModule.La * _numericsModule.Tscale;
            d2h_d2t = w * d2y_d2t * NumericsModule.La * _numericsModule.Tscale;
        }

        private void Soilsat(double t, out double y, out double dy_dt, out double d2y_d2t, int shg)
        {
            SplineHermite.SplineHermitVal(PropertiesModule.UnfrNdata, PropertiesModule.UnfrXdata[shg], PropertiesModule.UnfrC[shg], t, out y, out dy_dt, out d2y_d2t);
        }

        private void Tridiag(int m, int n0, int n1, IndexedArray<double> a, IndexedArray<double> b, IndexedArray<double> c, IndexedArray<double> d)
        {
            // Solution of a tridiagonal matrix

            for (var k = m + 1; k <= n1; ++k)
            {
                var xm = a[k] / b[k - 1];
                b[k] = b[k] - xm * c[k - 1];
                d[k] = d[k] - xm * d[k - 1];
            }

            d[n1] = d[n1] / b[n1];

            for (var i = m + 1; i <= n1; ++i)
            {
                var k = n1 + m - i;
                d[k] = (d[k] - c[k] * d[k + 1]) / b[k];
            }
        }

        public void UpdateProperties(double snow_height)
        {
            var dx = 0.0;
            var i = 0;
            
            // Provided the snow height, 
            // this code analyzes the grid, and inserts a new node representing the snow surface
            if (1.0e-5 < snow_height)
            {
                if (Math.Abs(_numericsModule.Xref[_numericsModule.Sx]) < snow_height)
                {
                    _numericsModule.X[_numericsModule.Sx] = -snow_height;
                    _numericsModule.Sn = _numericsModule.Sx;
                }
                else
                {
                    i = _numericsModule.Sx;
                    while (i < 0)
                    {
                        if ((Math.Abs(_numericsModule.Xref[i + 1]) < snow_height) && (snow_height <= Math.Abs(_numericsModule.Xref[i]) + 1.0e-10))
                        {
                            if ((Math.Abs(Math.Abs(_numericsModule.Xref[i + 1]) - snow_height) < Math.Abs(_numericsModule.Xref[i + 1] - _numericsModule.Xref[i]) * 0.5) && (i < -1))
                            {
                                _numericsModule.X[i + 1] = -snow_height;
                                _numericsModule.Sn = i + 1;
                            }
                            else
                            {
                                _numericsModule.X[i] = -snow_height;
                                _numericsModule.Sn = i;
                            }
                            i = 1;
                        }
                        else
                        {
                            i = i + 1;
                        }
                    }
                }
            }
            else
            {
                _numericsModule.Sn = 2;
            }

            for (var j = _numericsModule.Sx; j <= _numericsModule.Sn - 1; ++j)
            {
                _numericsModule.X[j] = _numericsModule.Xref[j];
            }

            for (var j = _numericsModule.Sn + 1; j <= 2; ++j)
            {
                _numericsModule.X[j] = _numericsModule.Xref[j];
            }

            for (var j = 2; j <= _numericsModule.Lx; ++j)
            {
                _numericsModule.X[j] = _numericsModule.Xref[j];
            }
            
            // Computes distances between each node in the grid in order to speed up consequent computations
            for (i = _numericsModule.Sn; i <= 1; ++i)
            {
                if (i == _numericsModule.Sn)
                {
                    dx = (_numericsModule.X[_numericsModule.Sn + 1] - _numericsModule.X[_numericsModule.Sn]) / 2.0;
                }
                else if (i == _numericsModule.Lx)
                {
                    dx = (_numericsModule.X[_numericsModule.Lx] - _numericsModule.X[_numericsModule.Lx - 1]) / 2.0;
                }
                else if ((i < _numericsModule.Lx) && (_numericsModule.Sx < i))
                {
                    dx = (_numericsModule.X[i] - _numericsModule.X[i - 1]) / 2.0 + (_numericsModule.X[i + 1] - _numericsModule.X[i]) / 2.0;
                }
                _numericsModule.Diag[i] = dx;
                if (i < _numericsModule.Lx) _numericsModule.Dx1[i] = 1.0 / (_numericsModule.X[i + 1] - _numericsModule.X[i]);
            }
            
            // Update the snow thermal properties
            for (i = _numericsModule.Sn; i <= 0; ++i)
            {
                // Heat capacity of ice
                _numericsModule.LayerCf[i] = _propertiesModule.SnowHeatc;
                // Thermal conductivity of snow
                _numericsModule.LayerLf[i] = _propertiesModule.SnowThcnd;
                // 
                // Water content=0,No phase change in snow
                _numericsModule.LayerW[i] = 0.0;
            }
            
            // Just for the sake of pretty looking graphics
            for (i = _numericsModule.Sx; i <= _numericsModule.Sn - 1; ++i)
            {
                // Heat capacity of ice
                _numericsModule.LayerCf[i] = -1.0;
                // Thermal conductivity of thawed soil
                _numericsModule.LayerLf[i] = -1.0;

                // Porosity=0,No phase change
                _numericsModule.LayerW[i] = 0.0;
            }
        }

        //public void InitializeSoil(out double soil_dt, out int gc, out double soil_time)        
        public void InitializeSoil()
        {
            // soilDt   Automatically updated timestep
            // gC       Automatically updated variable
            // soilTime Time in soil, for testing

            var dx = 0.0;

            // copy xref into x
            _numericsModule.X = _numericsModule.Xref.Copy();

            for (var i = _numericsModule.Sx; i <= _numericsModule.Lx; ++i)
            {
                if (i == _numericsModule.Sx)
                {
                    dx = (_numericsModule.X[_numericsModule.Sx + 1] - _numericsModule.X[_numericsModule.Sx]) / 2.0;
                }
                else if (i == _numericsModule.Lx)
                {
                    dx = (_numericsModule.X[_numericsModule.Lx] - _numericsModule.X[_numericsModule.Lx - 1]) / 2.0;
                }
                else if ((i < _numericsModule.Lx) && (_numericsModule.Sx < i))
                {
                    dx = (_numericsModule.X[i] - _numericsModule.X[i - 1]) / 2.0 + (_numericsModule.X[i + 1] - _numericsModule.X[i]) / 2.0;
                }
                _numericsModule.Diag[i] = dx;
                if (i < _numericsModule.Lx) _numericsModule.Dx1[i] = 1.0 / (_numericsModule.X[i + 1] - _numericsModule.X[i]);
            }

            // Update thermal properties at each soil layer

            for (var i = 1; i <= _numericsModule.Lx; ++i)
            {
                //_numericsModule.LayerW[i + 1] = _propertiesModule.SoilWaterContent[_numericsModule.LayerType[i + 1]];  //Volumetric water content
                _numericsModule.LayerLf[i + 1] = Math.Pow(PropertiesModule.SoilLm[_numericsModule.LayerType[i + 1]], 1.0 - _numericsModule.LayerW[i + 1]) * Math.Pow(NumericsModule.Ki, _numericsModule.LayerW[i + 1]);
                _numericsModule.LayerCf[i + 1] = PropertiesModule.SoilCm[_numericsModule.LayerType[i + 1]] * (1.0 - _numericsModule.LayerW[i + 1]) + NumericsModule.Ci * _numericsModule.LayerW[i + 1];
            }

            //soil_time = 0.0;
            //gc = 0;
            //soil_dt = NumericsModule.MaxDt;
            _numericsModule.Sn = 0;
        }
    }
}

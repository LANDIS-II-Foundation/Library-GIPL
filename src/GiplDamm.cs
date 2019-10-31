using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace Landis.Extension.GiplDamm
{
    public enum GiplDammLayerType
    {
        Snow = 0,
        LiveMoss = 1,
        DeadMoss = 2,
        Peat = 3,
        MineralSoil = 4
    }

    public partial class GiplDamm
    {
        #region fields
       
        private NumericsModule _numericsModule;
        private PropertiesModule _propertiesModule;

        private static Regex _whiteSpaceRegex = new Regex(@"\s+");

        private int _lx;
        private int _sx;

        private IndexedArray<double> _temperature;
        private IndexedArray<double> _liquidFraction;

        private string _wd;
        private int _instance;
        //const string _file = "ModeledGroundTemperatureGiplDamm";

        private List<double> _shawDepths;

        private static string ver = "1";


        // global static properties inputs
        private static HashSet<GiplDammLayerType> _enabledLayerTypes;
        private static double _wcrit;
        private static double _geothermalHeatFlux;
        //private static double _defaultSnowThermalConductivity;
        //private static double _defaultSnowVolumetricHeatCapacity;
        private static double[] _defaultSoilWaterContentByLayerType;

  
        #endregion

        private void Abort()
        {
            Console.WriteLine();
            Console.WriteLine(" ===> Simulation aborted; Press Enter to end");
            Console.ReadLine();
        }

        #region constructor
        public GiplDamm()
        {

        }
        #endregion

        #region properties

        public int Lx => _lx;
        public int Sx => _sx;

        public static bool HasGlobalSetup { get; private set; }

        public string Name { get; set; }

        #endregion

        #region methods

        public HashSet<GiplDammLayerType> EnabledLayerTypes => _enabledLayerTypes;

        public GiplDammLayerType[] GetNodeLayerTypes()
        {
            var layerTypes = new GiplDammLayerType[_lx];
            for (var i = 0; i < _lx; ++i)
                layerTypes[i] = (GiplDammLayerType)_numericsModule.LayerType[i + 2];

            return layerTypes;
        }

        public List<double> GetGiplDepths()
        {
            var depths = new List<double>(_numericsModule.Lx);

            for (var j = 1; j <= _numericsModule.Lx; ++j)
                depths.Add(_numericsModule.Xref[j]);

            return depths;
        }

        //public Dictionary<GiplDammLayerType, double> GetLayerTypePorosities()
        //{
        //    var porosities = new Dictionary<GiplDammLayerType, double>();
        //    foreach (GiplDammLayerType g in EnabledLayerTypes)
        //        porosities[g] = _propertiesModule.SoilWaterContent[(int)g];

        //    return porosities;
        //}

        public static bool GlobalInitialization(string wd)
        {
            HasGlobalSetup = true;

            List<string> lines;
            List<string> line;

            // **
            // read the 'properties.txt' file

            //if (!ReadTextFile(Path.Combine(wd, "propertiesGd.txt"), out lines))
            if (!ReadTextFile(Path.Combine(wd, "properties.txt"), out lines))
                return false;

            // initialize properties module global properties
            PropertiesModule.InitializeGlobalProperties();

            _enabledLayerTypes = new HashSet<GiplDammLayerType>();
            _defaultSoilWaterContentByLayerType = new double[Enum.GetNames(typeof(GiplDammLayerType)).Length];

            var rc = 0;

            ++rc;   // Geothermal heat flux header
            _geothermalHeatFlux = double.Parse(lines[rc++]);
            ++rc;   // W crit header
            _wcrit = double.Parse(lines[rc++]);
            //++rc;   // snow thermal conductivity header
            //_defaultSnowThermalConductivity = double.Parse(lines[rc++]);
            //++rc;   // snow volumetric heat capacity header
            //_defaultSnowVolumetricHeatCapacity = double.Parse(lines[rc++]);

            // skip blank rows
            while (string.IsNullOrEmpty(lines[rc]))
                ++rc;

            ++rc;   // layer properties table
            while (rc < lines.Count)
            {
                line = ParseLine(lines[rc++]);

                if (!Enum.TryParse(line[0], true, out GiplDammLayerType lType))
                {
                    Console.WriteLine($"Unrecognized GiplDamm Layer Type '{line[0]}' in properties file");
                    return false;
                }

                _enabledLayerTypes.Add(lType);  // keep track of the layer types

                var i = (int)lType;
                _defaultSoilWaterContentByLayerType[i] = double.Parse(line[1]);
                PropertiesModule.SoilLm[i] = double.Parse(line[2]);
                PropertiesModule.SoilCm[i] = double.Parse(line[3]);

                // last column refers to another file
                List<string> flines;
                if (!ReadTextFile(Path.Combine(wd, line[4]), out flines))
                    return false;

                var frc = 0;   // file row counter
                for (var j = 1; j <= PropertiesModule.UnfrNdata; ++j)
                {
                    var fline = _whiteSpaceRegex.Split(flines[frc++]);
                    PropertiesModule.UnfrXdata[i][j] = double.Parse(fline[0]);
                    PropertiesModule.UnfrFdata[i][j] = double.Parse(fline[1]);
                    PropertiesModule.UnfrDfdata[i][j] = double.Parse(fline[2]);
                }

                // create Hermite polynomial
                SplineHermite.SplineHermiteSet(PropertiesModule.UnfrNdata, PropertiesModule.UnfrXdata[i], PropertiesModule.UnfrFdata[i], PropertiesModule.UnfrDfdata[i], PropertiesModule.UnfrC[i]);
            }

            return true;
        }

        public bool Initialize(string wd, string name, List<double> shawDepths)
        {
            _wd = wd;
            Name = name;

            //_instance = instance;
            _shawDepths = shawDepths;

            List<string> lines;
            List<string> line;

            _propertiesModule = new PropertiesModule();


            // **
            // read the 'initial.txt' file

            if (!ReadTextFile(Path.Combine(_wd, $"initial.txt"), out lines))
            //if (!ReadTextFile(Path.Combine(_wd, $"initialGd{_instance}.txt"), out lines))
                return false;

            var rc = 0;

            line = ParseLine(lines[rc++]);
            _lx = int.Parse(line[0]);
            _sx = int.Parse(line[1]);

            _temperature = new IndexedArray<double>(_sx, _lx + 1);
            _liquidFraction = new IndexedArray<double>(_sx, _lx + 1);

            _numericsModule = new NumericsModule(_sx, _lx);

            for (var i = _sx; i <= 0; ++i)
            {
                line = ParseLine(lines[rc++]);
                //_numericsModule.LayerType[i] = int.Parse(line[1]);
                _numericsModule.Xref[i] = double.Parse(line[2]);
                _temperature[i] = double.Parse(line[3]);
            }


            for (var i = 1; i <= _lx; ++i)
            {
                line = ParseLine(lines[rc++]);

                if (!Enum.TryParse(line[1], true, out GiplDammLayerType lType) || !_enabledLayerTypes.Contains(lType))
                {
                    Console.WriteLine($"Unrecognized GiplDamm Layer Type '{line[1]}' or the Layer Type not enabled in the properties file, at soil layer {i}");
                    return false;
                }

                _numericsModule.LayerType[i + 1] = (int)lType;
                _numericsModule.Xref[i] = double.Parse(line[2]);
                _temperature[i + 1] = double.Parse(line[3]);

                // set the default water content fraction profile based on the default layer type water content
                _numericsModule.DefaultWaterContentProfile[i + 1] = _defaultSoilWaterContentByLayerType[(int)lType];
            }


            // initialize liquid fraction

            for (var i = 2; i <= _numericsModule.Lx + 1; ++i)
            {
                double y, dummy1, dummy2;
                Soilsat(_temperature[i], out y, out dummy1, out dummy2, _numericsModule.LayerType[i]);
                _liquidFraction[i] = y;
            }

            for (var i = _numericsModule.Sn; i <= 0; ++i)
            {
                _liquidFraction[i] = 0.0;
            }


            //var writer = new StreamWriter(Path.Combine(_wd, $"{_file}{_instance}.txt"));
            //var writer = new StreamWriter(Path.Combine(_wd, $"{_file}.txt"));

            //var time = 0.0;
            //writer.Write($"{time,15:F8}");
            //for (var j = 1; j <= _numericsModule.Lx; ++j)
            //    writer.Write($"{_numericsModule.Xref[j],15:F8}");
            //writer.WriteLine();

            //writer.Close();

            var writer = new StreamWriter(Path.Combine(_wd, $"{Name}_Log.txt"));

            var time = 0.0;
            writer.Write($"{time,15:F8}");
            foreach (var d in shawDepths)
                writer.Write($"{d,15:F8}");
            writer.WriteLine();

            writer.Close();

            return true;
        }

        public GiplDammResults CalculateSoilTemperature(double[] dailyAirTemp, double[] dailySnowThickness, double[] shawWaterContentFraction,
            double snowThermalConductivity, double snowVolumetricHeatCapacity)
        {
            // units:
            //  dailySnowThickness          m
            //  shawDepths                  m
            //  snowThermalConductivity     W/m/K
            //  snowVolumetricHeatCapacity  J/m3/K

            var results = new GiplDammResults();

            // overwrite water content if input is not null.  
            //  otherwise, use the default water content profile based on the properties file.  (This should only happen for the first call to Gipl.)
            if (shawWaterContentFraction != null)
            {
                // the Shaw water content profile is at the shaw depths.
                //  map these onto the Gipl depth points:  
                //   if Gipl depth < the first Shaw depth, use the value at the Shaw depth,
                //   else if Gipl depth > the last Shaw depth, use the value at the last Shaw depth.
                //   otherwise, find the Shaw depths that bracket the Gipl depth and interpolate the values.

                var j = 0;
                for (var i = 1; i <= _numericsModule.Lx; ++i)
                {
                    // find the first water content depth index j such that waterContentDepth[j] > Xref[i];
                    while (j < _shawDepths.Count && _shawDepths[j] < _numericsModule.Xref[i])
                        ++j;

                    if (j == 0)
                        _numericsModule.LayerW[i + 1] = shawWaterContentFraction.First();
                    else if (j == _shawDepths.Count)
                        _numericsModule.LayerW[i + 1] = shawWaterContentFraction.Last();
                    else
                        _numericsModule.LayerW[i + 1] = shawWaterContentFraction[j - 1] + (shawWaterContentFraction[j] - shawWaterContentFraction[j - 1]) * (_numericsModule.Xref[i] - _shawDepths[j - 1]) / (_shawDepths[j] - _shawDepths[j - 1]);
                }
            }
            else
                for (var i = 1; i <= _numericsModule.Lx; ++i)
                {
                    _numericsModule.LayerW[i + 1] = _numericsModule.DefaultWaterContentProfile[i + 1];
                }

            _propertiesModule.SnowThcnd = snowThermalConductivity;

            // GIPL expects snow volumetric heat capacity to be in MJ/m3/K, but so I need to convert the passed value, which is in J/m3/K
            _propertiesModule.SnowHeatc = snowVolumetricHeatCapacity * 1.0e-6;

            var npoints = dailyAirTemp.Length;

            InitializeSoil();

            var snowHeight = dailySnowThickness[0];

            //var writer = new StreamWriter(Path.Combine(_wd, $"{_file}{_instance}.txt"), true);
            //var writer = new StreamWriter(Path.Combine(_wd, $"{_file}.txt"), true);
            var writer = new StreamWriter(Path.Combine(_wd, $"{Name}_Log.txt"), true);

            var time = 0.0;         // automatically update time
            var dtime = NumericsModule.MaxDt;   // automatically updated timestep
            var gc = 0;                                  // automatically updated counter of converged iterations

            for (var i = 0; i < npoints; ++i)
            {
                // John McNabb: timestep is daily, so timeEnd is always equal to i.
                //var timeEnd = _bcModule.Timepoints[i][1];
                double timeEnd = i;
                snowHeight = dailySnowThickness[i];

                var surfaceTemp = SurfaceTemperature(timeEnd, dailyAirTemp);
                //Console.WriteLine($"t={timeEnd,6:F1}; AirT={surfaceTemp,5:F1}; SnowHGT={snowHeight,5:F1}; Ground Temp={_temperature[71],6:F2} @ Z={_numericsModule.X[71],4:F2}");

                UpdateProperties(snowHeight);

                Stemperature(_temperature, _liquidFraction, ref time, timeEnd, ref dtime, ref gc, dailyAirTemp);

                //writer.Write($"{time,15:F8}");
                //for (var j = 2; j <= _numericsModule.Lx + 1; ++j)
                //    writer.Write($"{_temperature[j],15:F8}");
                //writer.WriteLine();

                // save the temperature profile
                var soilTempProfile = new double[_numericsModule.Lx];
                for (var jj = 2; jj <= _numericsModule.Lx + 1; ++jj)
                    soilTempProfile[jj - 2] = _temperature[jj];

                results.DailySoilTemperatureProfiles.Add(soilTempProfile);

                if (_shawDepths != null)
                {
                    // interpolate temperature results onto Shaw depths for results.
                    //  use the reverse logic from what was used above to map the Shaw water content profile
                    //  onto the Gipl depths.

                    var shawTemps = new double[_shawDepths.Count];

                    var j = 1;
                    for (var k = 0; k < _shawDepths.Count; ++k)
                    {
                        // find the first index j such that Xref[j] > shawDepths[k].
                        // Xref[j] corresponds to _temperature[j + 1].
                        while (j < _numericsModule.Lx + 1 && _numericsModule.Xref[j] < _shawDepths[k])
                            ++j;

                        if (j == 1)
                            shawTemps[k] = _temperature[2];
                        else if (j == _numericsModule.Lx + 1)
                            shawTemps[k] = _temperature[_numericsModule.Lx + 1];
                        else
                            shawTemps[k] = _temperature[j] + (_temperature[j + 1] - _temperature[j]) * (_shawDepths[k] - _numericsModule.Xref[j - 1]) / (_numericsModule.Xref[j] - _numericsModule.Xref[j - 1]);

                    }

                    results.DailySoilTemperatureProfilesAtShawDepths.Add(shawTemps);

                    writer.Write($"{time,15:F8}");
                    foreach (var t in shawTemps)
                        writer.Write($"{t,15:F8}");
                    writer.WriteLine();

                }
            }

            writer.Close();

            results.MakeProfileAveragesOverDays();
            results.Success = true;
            return results;
        }

        #endregion

        #region private methods

        private double SurfaceTemperature(double time, double[] airTemp)
        {
            double surface_temperature;

            // John McNabb: assume timesteps are always daily so dtimepoints = 1.0;
            //var dtimepoints = _bcModule.Timepoints[1][1] - _bcModule.Timepoints[0][1];
            var dtimepoints = 1.0;
            var i = (int)(time / dtimepoints);   // FORTRAN code uses dint() which truncates
            var di = time - i * dtimepoints;

            if (i < airTemp.Length - 1)
            {
                surface_temperature = airTemp[i] + di * (airTemp[i + 1] - airTemp[i]) / dtimepoints;
            }
            else
            {
                surface_temperature = airTemp.Last();
            }

            return surface_temperature;
        }

        private static bool ReadTextFile(string filePath, out List<string> lines)
        {
            lines = null;
            if (!File.Exists(filePath))
            {
                Console.WriteLine($"Cannot find file: '{filePath}'.");
                return false;
            }

            lines = File.ReadAllLines(filePath).Select(x => x.Trim()).ToList();

            // remove blank lines at the end of the file
            while (string.IsNullOrEmpty(lines.Last()))
                lines.RemoveAt(lines.Count - 1);

            if (!lines.Any())
            {
                Console.WriteLine($"File: '{filePath}' is empty.");
                return false;
            }

            return true;
        }

        private static List<string> ParseLine(string s)
        {
            return _whiteSpaceRegex.Split(s.Trim()).ToList();
        }

        #endregion
    }
}

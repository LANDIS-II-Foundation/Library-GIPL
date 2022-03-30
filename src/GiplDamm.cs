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

        private static string _wd;

        private int _maxCommonNode;

        // global static properties inputs
        private static double _wcrit;
        private static double _geothermalHeatFlux;
        private static List<string> _layerTypes;

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

        public static double SnowC1 { get; set; }
        
        #endregion

        #region methods

        public List<double> GetGiplDepths()
        {
            var depths = new List<double>(_numericsModule.Lx);

            for (var j = 1; j <= _numericsModule.Lx; ++j)
                depths.Add(_numericsModule.Xref[j]);

            return depths;
        }

        public static bool GlobalInitialization(string propertiesFilePath, out string errorMessage)
        {
            errorMessage = string.Empty;
            _wd = Path.GetDirectoryName(propertiesFilePath);

            HasGlobalSetup = true;

            List<string> lines;
            List<string> line;

            // **
            // read the 'properties.txt' file

            if (!ReadTextFile(propertiesFilePath, out lines, out errorMessage))
                return false;

            var rc = -1;

            _geothermalHeatFlux = _wcrit = double.NaN;
            SnowC1 = 8.33333E-08;   // default value equal to 0.03e-2 / 3600

            var headersFound = false;

            while (++rc < lines.Count)
            {
                if (lines[rc].StartsWith("LayerType", StringComparison.OrdinalIgnoreCase))
                {
                    headersFound = true;
                    break;
                }

                if (lines[rc].StartsWith("GeothermalHeatFlux", StringComparison.OrdinalIgnoreCase))
                {
                    line = ParseLine(lines[rc]);
                    if (!double.TryParse(line[1], out _geothermalHeatFlux))
                    {
                        errorMessage = $"Cannot parse '{line[1]}' as numeric for GeothermalHeatFlux";
                        return false;
                    }
                }

                if (lines[rc].StartsWith("WCritial", StringComparison.OrdinalIgnoreCase))
                {
                    line = ParseLine(lines[rc]);
                    if (!double.TryParse(line[1], out _wcrit))
                    {
                        errorMessage = $"Cannot parse '{line[1]}' as numeric for WCritial";
                        return false;
                    }
                }

                if (lines[rc].StartsWith("SnowC1", StringComparison.OrdinalIgnoreCase))
                {
                    line = ParseLine(lines[rc]);
                    if (!double.TryParse(line[1], out var c1))
                    {
                        errorMessage = $"Cannot parse '{line[1]}' as numeric for SnowC1";
                        return false;
                    }
                    SnowC1 = c1;
                }
            }

            if (!headersFound)
            {
                errorMessage = $"Cannot find properties table containing 'LayerType', etc.";
                return false;
            }

            if (double.IsNaN(_geothermalHeatFlux) || double.IsNaN(_wcrit))
            {
                errorMessage = $"Cannot find 'GeothermalHeatFlux' and/or 'WCritial'";
                return false;
            }

            // parse layer properties table
            ++rc;   // skip headers
            var layerRows = new List<string>();

            while (rc < lines.Count && !string.IsNullOrEmpty(lines[rc]))
            {
                layerRows.Add(lines[rc]);
                ++rc;
            }


            // initialize properties module global properties
            PropertiesModule.InitializeGlobalProperties(layerRows.Count);

            _layerTypes = new List<string>();
            var i = 1;
            foreach (var layerRow in layerRows)
            {
                line = ParseLine(layerRow);
                _layerTypes.Add(line[0]);

                PropertiesModule.SoilLm[i] = double.Parse(line[1]);
                PropertiesModule.SoilCm[i] = double.Parse(line[2]);

                // last column refers to another file
                List<string> flines;
                if (!ReadTextFile(Path.Combine(_wd, line[3]), out flines, out errorMessage))
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

                ++i;
            }

            return true;
        }

        public bool Initialize(string name, Dictionary<string, string> thuFileData, List<double> soilThetaSat, out string errorMessage)
        {
            errorMessage = string.Empty;

            Name = name;
            _propertiesModule = new PropertiesModule();

            // get the layers to be processed

            var snowLayerInput = new GiplLayerInput();

            if (!TryParseInput("SnowNodes", thuFileData["SnowNodes"], out snowLayerInput.Nodes, out errorMessage, 0))
                return false;

            if (!TryParseInput("MaxSnowThickness", thuFileData["MaxSnowThickness"], out snowLayerInput.MaxDepth, out errorMessage, 0.0))
                return false;

            if (!TryParseInput("InitSnowTemperature", thuFileData["InitSnowTemperature"], out snowLayerInput.InitTemperature, out errorMessage, AbsZero))
                return false;

            // read common layers 
            var soilLayerInputs = new List<GiplLayerInput>();
            var jj = 1;
            while (true)
            {
                var key = $"GiplSoilType{jj}";
                if (!thuFileData.TryGetValue(key, out var soilType) || string.IsNullOrEmpty(soilType))
                    break;

                if (!_layerTypes.Contains(soilType, StringComparer.OrdinalIgnoreCase))
                {
                    errorMessage = $"Unrecognized {key} '{soilType}'. Available types are {string.Join(",", _layerTypes)}";
                    return false;
                }

                var layer = new GiplLayerInput { LayerIndex = _layerTypes.FindIndex(x => x.Equals(soilType, StringComparison.OrdinalIgnoreCase)) + 1 };

                if (!TryParseInput($"Nodes{jj}", thuFileData[$"Nodes{jj}"], out layer.Nodes, out errorMessage, 0))
                    return false;

                if (!TryParseInput($"MaxDepth{jj}", thuFileData[$"MaxDepth{jj}"], out layer.MaxDepth, out errorMessage, 0.0))
                    return false;

                if (!TryParseInput($"InitTemperature{jj}", thuFileData[$"InitTemperature{jj}"], out layer.InitTemperature, out errorMessage, AbsZero))
                    return false;

                if (!TryParseInput($"InitWaterContent{jj}", thuFileData[$"InitWaterContent{jj}"], out layer.InitWaterContent, out errorMessage, 0.0))
                    return false;

                if (soilLayerInputs.Any() && layer.MaxDepth <= soilLayerInputs.Last().MaxDepth)
                {
                    errorMessage = $"MaxDepth for layer {jj} ({layer.MaxDepth}) is <= the MaxDepth ({soilLayerInputs.Last().MaxDepth}) of the previous layer.  MaxDepths must increase for each Gipl layer.";
                    return false;
                }

                soilLayerInputs.Add(layer);
                ++jj;
            }

            _maxCommonNode = soilLayerInputs.Sum(x => x.Nodes) + 1;

            // read bottom layer (optional)
            if (!string.IsNullOrEmpty(thuFileData["BottomGiplSoilType"]))
            {
                var soilType = thuFileData["BottomGiplSoilType"];
                if (!_layerTypes.Contains(soilType, StringComparer.OrdinalIgnoreCase))
                {
                    errorMessage = $"Unrecognized BottomGiplSoilType '{soilType}'. Available types are {string.Join(",", _layerTypes)}";
                    return false;
                }

                var layer = new GiplLayerInput { LayerIndex = _layerTypes.FindIndex(x => x.Equals(soilType, StringComparison.OrdinalIgnoreCase)) + 1};

                if (!TryParseInput("BottomGiplNodes", thuFileData["BottomGiplNodes"], out layer.Nodes, out errorMessage, 0))
                    return false;

                if (!TryParseInput("BottomGiplMaxDepth", thuFileData["BottomGiplMaxDepth"], out layer.MaxDepth, out errorMessage, 0.0))
                    return false;

                if (!TryParseInput("BottomGiplInitTemperature", thuFileData["BottomGiplInitTemperature"], out layer.InitTemperature, out errorMessage, AbsZero))
                    return false;

                if (!TryParseInput("BottomGiplWaterContent", thuFileData["BottomGiplWaterContent"], out layer.InitWaterContent, out errorMessage, 0.0))
                    return false;

                if (soilLayerInputs.Any() && layer.MaxDepth <= soilLayerInputs.Last().MaxDepth)
                {
                    errorMessage = $"BottomGiplMaxDepth ({layer.MaxDepth}) is <= the MaxDepth ({soilLayerInputs.Last().MaxDepth}) of the previous layer.  MaxDepths must increase for each Gipl layer.";
                    return false;
                }

                soilLayerInputs.Add(layer);
            }

            // set up layers
            _sx = -snowLayerInput.Nodes;
            _lx = soilLayerInputs.Sum(x => x.Nodes) + 1;

            _temperature = new IndexedArray<double>(_sx, _lx + 1);
            _liquidFraction = new IndexedArray<double>(_sx, _lx + 1);
            _numericsModule = new NumericsModule(_sx, _lx);

            // snow
            var step = snowLayerInput.MaxDepth / snowLayerInput.Nodes;

            for (var i = _sx; i <= 0; ++i)
            {
                _numericsModule.Xref[i] = step * i;
                _temperature[i] = snowLayerInput.InitTemperature;
            }

            var previousMaxDepth = 0.0;
            var j = 1;
            foreach (var layer in soilLayerInputs)
            {
                var maxDepth = layer.MaxDepth;
                step = (maxDepth - previousMaxDepth) / layer.Nodes;
                for (var i = layer == soilLayerInputs.First() ? 0 : 1; i <= layer.Nodes; ++i)
                {
                    _numericsModule.LayerType[j + 1] = layer.LayerIndex;
                    _numericsModule.Xref[j] = previousMaxDepth + step * i;
                    _temperature[j + 1] = layer.InitTemperature;

                    // set the default water content fraction profile to the initial water content
                    _numericsModule.LayerW[j + 1] = layer.InitWaterContent;

                    ++j;
                }

                previousMaxDepth = maxDepth;
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

            var logPath = Path.Combine(_wd, "Outputs");
            if (!Directory.Exists(logPath))
                Directory.CreateDirectory(logPath);

            var writer = new StreamWriter(Path.Combine(logPath, $"Log_{Name}.txt"));

            var year = 0;
            var time = 0.0;
            writer.Write($"{year,6}{time,15:F8}");
            //for (var i = 0; i < _maxCommonNode; ++i)
            for (var i = 0; i < _numericsModule.Lx; ++i)
                writer.Write($"{_numericsModule.Xref[i + 1],15:F8}");

            writer.WriteLine();

            writer.Close();

            var snowWriter = new StreamWriter(Path.Combine(logPath, $"Log_{Name}_Snow.csv"));
            snowWriter.WriteLine("Year,Month,Day,SnowThickness[m],SnowThermCond[W/m/K],SnowVolHeatCap[J/m3/K]");
            snowWriter.Close();

            // set LayerP to the soil theta sat from Shaw
            for (var i = 1; i <= _numericsModule.Lx; ++i)
            {
                _numericsModule.LayerP[i + 1] = (i - 1 >= soilThetaSat.Count) ? soilThetaSat.Last() : soilThetaSat[i - 1];
            }

            return true;
        }

        private class GiplLayerInput
        {
            public int LayerIndex;
            public int Nodes;
            public double MaxDepth;
            public double InitTemperature;
            public double InitWaterContent;
        }

        public GiplDammResults CalculateSoilTemperature(int year, int month, double[] dailyAirTemp, double[] dailySnowThickness, double[] dailySnowThermalConductivities, double[] dailySnowVolumetricHeatCapacities, double[] shawTotalWaterContentFraction)
        {
            // units:
            //  dailySnowThickness                 m
            //  dailySnowThermalConductivities     W/m/K
            //  dailySnowVolumetricHeatCapacities  J/m3/K

            var snowWriter = new StreamWriter(Path.Combine(_wd, "Outputs", $"Log_{Name}_Snow.csv"), true);
            for (var day = 0; day < dailySnowThickness.Length; ++day)
                snowWriter.WriteLine($"{year},{month},{day},{dailySnowThickness[day]},{dailySnowThermalConductivities[day]},{dailySnowVolumetricHeatCapacities[day]}");
            snowWriter.Close();

            // calculate snow thermal conductivity and volumetric heat capacity from the non-zero daily values
            var nonZeroTk = dailySnowThermalConductivities.Where(x => x > 0.0).ToList();
            var nonZeroHc = dailySnowVolumetricHeatCapacities.Where(x => x > 0.0).ToList();
            var snowThermalConductivity = nonZeroTk.Any() ? nonZeroTk.Average() : 0.0;
            var snowVolumetricHeatCapacity = nonZeroHc.Any() ? nonZeroHc.Average() : 0.0;

            var results = new GiplDammResults();

            // overwrite water content if input is not null.  
            //  otherwise, use the default water content profile based on the properties file.  (This should only happen for the first call to Gipl.)
            if (shawTotalWaterContentFraction != null)
            {
                for (var i = 0; i < _maxCommonNode; ++i)
                {
                    // use the minimum of the Shaw total water content and LayerP
                    _numericsModule.LayerW[i + 2] = Math.Min(shawTotalWaterContentFraction[i], _numericsModule.LayerP[i + 2]);
                }
            }

            _propertiesModule.SnowThcnd = snowThermalConductivity;

            // GIPL expects snow volumetric heat capacity to be in MJ/m3/K, but so I need to convert the passed value, which is in J/m3/K
            _propertiesModule.SnowHeatc = snowVolumetricHeatCapacity * 1.0e-6;

            var npoints = dailyAirTemp.Length;

            InitializeSoil();

            var snowHeight = dailySnowThickness[0];

            var writer = new StreamWriter(Path.Combine(_wd, "Outputs", $"Log_{Name}.txt"), true);

            var time = 0.0;         // automatically update time
            var dtime = NumericsModule.MaxDt;   // automatically updated timestep
            var gc = 0;                                  // automatically updated counter of converged iterations

            for (var i = 0; i < npoints; ++i)
            {
                // John McNabb: timestep is daily, so timeEnd is always equal to i.
                double timeEnd = i;
                snowHeight = dailySnowThickness[i];

                var surfaceTemp = SurfaceTemperature(timeEnd, dailyAirTemp);
                //Console.WriteLine($"t={timeEnd,6:F1}; AirT={surfaceTemp,5:F1}; SnowHGT={snowHeight,5:F1}; Ground Temp={_temperature[71],6:F2} @ Z={_numericsModule.X[71],4:F2}");

                UpdateProperties(snowHeight);

                Stemperature(_temperature, _liquidFraction, ref time, timeEnd, ref dtime, ref gc, dailyAirTemp);

                // save the temperature profile
                var soilTempProfile = new double[_numericsModule.Lx];
                for (var jj = 2; jj <= _numericsModule.Lx + 1; ++jj)
                    soilTempProfile[jj - 2] = _temperature[jj];

                results.DailySoilTemperatureProfiles.Add(soilTempProfile);

                var shawTemps = new double[_maxCommonNode];

                for (var k = 0; k < _maxCommonNode; ++k)
                    shawTemps[k] = _temperature[2 + k];

                results.DailySoilTemperatureProfilesAtShawDepths.Add(shawTemps);

                writer.Write($"{year,6}{time,15:F8}");
                //foreach (var t in shawTemps)
                for (var k = 0; k < _numericsModule.Lx; ++k)
                    writer.Write($"{_temperature[2 + k],15:F8}");
                writer.WriteLine();
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

        private static bool ReadTextFile(string filePath, out List<string> lines, out string errorMessage)
        {
            errorMessage = string.Empty;
            lines = null;
            if (!File.Exists(filePath))
            {
                errorMessage = $"Cannot find file: '{filePath}'.";
                return false;
            }

            lines = File.ReadAllLines(filePath).Select(x => x.Trim()).ToList();

            // remove blank lines at the end of the file
            while (string.IsNullOrEmpty(lines.Last()))
                lines.RemoveAt(lines.Count - 1);

            if (!lines.Any())
            {
                errorMessage = $"File: '{filePath}' is empty.";
                return false;
            }

            return true;
        }

        private static List<string> ParseLine(string s)
        {
            return _whiteSpaceRegex.Split(s.Trim()).ToList();
        }

        private static double AbsZero = -273.15;

        public static bool TryParseInput<T>(string name, string input, out T value, out string errorMessage, double lowerRange = double.NegativeInfinity, double upperRange = double.PositiveInfinity, bool lowerInclusive = true, bool upperInclusive = true)
        {
            value = default(T);

            if (string.IsNullOrEmpty(input))
            {
                errorMessage = $"{name} : input value is empty";
                return false;
            }

            errorMessage = string.Empty;
            if (typeof(T) == typeof(string))
            {
                value = (T)Convert.ChangeType(input, typeof(T));
                return true;
            }

            if (typeof(T) == typeof(double))
            {
                double t;
                if (!double.TryParse(input, out t))
                {
                    errorMessage = $"{name} : cannot parse '{input}' as double";
                    return false;
                }

                if (!CheckRange(t, lowerRange, upperRange, lowerInclusive, upperInclusive, out errorMessage))
                    return false;

                value = (T)Convert.ChangeType(t, typeof(T));
                return true;
            }

            if (typeof(T) == typeof(int))
            {
                int t;
                if (!int.TryParse(input, out t))
                {
                    errorMessage = $"{name} : cannot parse '{input}' as int";
                    return false;
                }

                if (!CheckRange(t, lowerRange, upperRange, lowerInclusive, upperInclusive, out errorMessage))
                    return false;

                value = (T)Convert.ChangeType(t, typeof(T));
                return true;
            }

            if (typeof(T) == typeof(bool))
            {
                bool t;
                if (!bool.TryParse(input, out t))
                {
                    errorMessage = $"{name} : cannot parse '{input}' as bool";
                    return false;
                }

                value = (T)Convert.ChangeType(t, typeof(T));
                return true;
            }

            errorMessage = $"{name} : unrecognized Type '{typeof(T)}' requested";

            return false;
        }

        public static bool CheckRange(double t, double lowerRange, double upperRange, bool lowerInclusive, bool upperInclusive, out string errorMessage)
        {
            errorMessage = string.Empty;

            var inRange = t > lowerRange && t < upperRange;

            inRange |= lowerInclusive && t == lowerRange;
            inRange |= upperInclusive && t == upperRange;

            if (!inRange)
            {
                errorMessage = $"Value '{t}' is out of the range: {(lowerInclusive && !double.IsNegativeInfinity(lowerRange) ? "[" : "(")}{lowerRange}, {upperRange}{(upperInclusive && !double.IsPositiveInfinity(upperRange) ? "]" : ")")}";
                return false;
            }

            return true;
        }

        #endregion
    }
}

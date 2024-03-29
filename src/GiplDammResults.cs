﻿using System;
using System.Collections.Generic;
using System.Text;
using System.Linq;

namespace Landis.Extension.GiplDamm
{
    public class GiplDammResults
    {
        public GiplDammResults()
        {
            DailySoilTemperatureProfilesAtShawDepths = new List<double[]>();
            DailySoilTemperatureProfiles = new List<double[]>();
        }

        public bool Success { get; set; }

        public List<double[]> DailySoilTemperatureProfilesAtShawDepths { get; }
        public List<double[]> DailySoilTemperatureProfiles { get; }

        public List<double> AverageSoilTemperatureProfileAtShawDepths { get; private set; }     // average across the month
        public List<double> AverageSoilTemperatureProfile { get; private set; }                 // average across the month

        public void MakeProfileAveragesOverDays()
        {
            AverageSoilTemperatureProfileAtShawDepths = AverageProfileOverDays(DailySoilTemperatureProfilesAtShawDepths);
            AverageSoilTemperatureProfile = AverageProfileOverDays(DailySoilTemperatureProfiles);
        }

        private List<double> AverageProfileOverDays(List<double[]> dailyProfiles)
        {
            var days = dailyProfiles.Count;
            var depths = dailyProfiles.First().Length;
            var averageProfile = new List<double>(depths);

            for (var j = 0; j < depths; ++j)
                averageProfile.Add(Enumerable.Range(0, days).Average(i => dailyProfiles[i][j]));

            return averageProfile;
        }
    }
}

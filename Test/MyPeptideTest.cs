﻿using Chemistry;
using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.Indexing;
using EngineLayer.ModernSearch;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public static class MyPeptideTest
    {
        [Test]
        public static void TestIdenticalPeaks()
        {
            IDictionary<int, List<Modification>> mods = new Dictionary<int, List<Modification>>();
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif);
            mods.Add(1, new List<Modification> { new ModificationWithMass("Hehe", null, motif, TerminusLocalization.NProt, 18.010565, null, null, null, null) });
            var prot = new Protein("MMMM", null, null, null, mods);
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 1);
            var ye = prot.Digest(digestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass>()).First();

            var massArray = ye.CompactPeptide(TerminusType.None).ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { ProductType.B, ProductType.Y });
            Array.Sort(massArray);
            double[] intensities = new double[] { 1, 1, 1, 1 };
            double[] mz = new double[] { massArray[0].ToMz(1), massArray[2].ToMz(1), massArray[4].ToMz(1), 10000 };
            MzSpectrum massSpectrum = new MzSpectrum(mz, intensities, false);
            MsDataScan scan = new MsDataScan(massSpectrum, 1, 1, true, Polarity.Positive, 1, new MzRange(300, 2000), "", MZAnalyzerType.Unknown, massSpectrum.SumOfAllY, null, null, "scan=1", 0, null, null, 0, null, DissociationType.Unknown, 1, null);

            PeptideSpectralMatch[] globalPsms = new PeptideSpectralMatch[1];
            Ms2ScanWithSpecificMass[] arrayOfSortedMS2Scans = { new Ms2ScanWithSpecificMass(scan, 0, 0, null) };
            CommonParameters CommonParameters = new CommonParameters(
                productMassTolerance: new PpmTolerance(5),
                scoreCutoff: 1,
                digestionParams: new DigestionParams(
                    maxMissedCleavages: 0,
                    minPeptideLength: 1,
                    initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain));
            ClassicSearchEngine cse = new ClassicSearchEngine(globalPsms, arrayOfSortedMS2Scans, new List<ModificationWithMass>(), new List<ModificationWithMass>(), new List<Protein> { prot }, new List<ProductType> { ProductType.B, ProductType.Y }, new OpenSearchMode(), CommonParameters, new List<string>());

            cse.Run();
            Assert.AreEqual(globalPsms[0].MatchedFragmentIons.Count, 3);
        }

        [Test]
        public static void TestLastPeaks()
        {
            IDictionary<int, List<Modification>> mods = new Dictionary<int, List<Modification>>();
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif);
            var prot = new Protein("MMMM", null, null, null, mods);
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 1);
            PeptideWithSetModifications thePep = prot.Digest(digestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass>()).First();
            PepsWSMods pep = new PepsWSMods(thePep.Protein, thePep.DigestionParams, thePep.OneBasedStartResidueInProtein, thePep.OneBasedEndResidueInProtein, thePep.PeptideDescription, thePep.MissedCleavages, thePep.AllModsOneIsNterminus, thePep.NumFixedMods);

            var massArray = pep.CompactPeptide(TerminusType.None).ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { ProductType.B, ProductType.Y });
            Array.Sort(massArray);
            double[] intensities = new double[] { 1, 1, 1 };
            double[] mz = new double[] { 1, 2, massArray[4].ToMz(1) };
            MzSpectrum massSpectrum = new MzSpectrum(mz, intensities, false);
            MsDataScan scan = new MsDataScan(massSpectrum, 1, 1, true, Polarity.Positive, 1, new MzRange(300, 2000), "", MZAnalyzerType.Unknown, massSpectrum.SumOfAllY, null, null, "scan=1", 0, null, null, 0, null, DissociationType.Unknown, 1, null);

            PeptideSpectralMatch[] globalPsms = new PeptideSpectralMatch[1];
            Ms2ScanWithSpecificMass[] arrayOfSortedMS2Scans = { new Ms2ScanWithSpecificMass(scan, 0, 0, null) };
            CommonParameters CommonParameters = new CommonParameters(
                scoreCutoff: 1,
                productMassTolerance: new PpmTolerance(5),
                digestionParams: new DigestionParams(
                    maxMissedCleavages: 0,
                    minPeptideLength: 1,
                    maxModificationIsoforms: int.MaxValue,
                    initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain));
            ClassicSearchEngine cse = new ClassicSearchEngine(globalPsms, arrayOfSortedMS2Scans, new List<ModificationWithMass>(), new List<ModificationWithMass>(), new List<Protein> { prot }, new List<ProductType> { ProductType.B, ProductType.Y }, new OpenSearchMode(), CommonParameters, new List<string>());

            cse.Run();
            Assert.Less(globalPsms[0].Score, 2);
            Assert.Greater(globalPsms[0].Score, 1);
        }

        [Test]
        public static void TestVeryCloseExperimentalsClassic()
        {
            IDictionary<int, List<Modification>> mods = new Dictionary<int, List<Modification>>();
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif);
            var prot = new Protein("MMMM", null, null, null, mods);
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 1);
            var thePep = prot.Digest(digestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass>()).First();

            var massArray = thePep.CompactPeptide(TerminusType.None).ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { ProductType.B, ProductType.Y });
            Array.Sort(massArray);
            double[] intensities = new double[] { 1, 1, 1, 1 };
            double[] mz = new double[] { 1, 2, massArray[4].ToMz(1), massArray[4].ToMz(1) + 1e-9 };
            MzSpectrum massSpectrum = new MzSpectrum(mz, intensities, false);
            MsDataScan scan = new MsDataScan(massSpectrum, 1, 1, true, Polarity.Positive, 1, new MzRange(300, 2000), "", MZAnalyzerType.Unknown, massSpectrum.SumOfAllY, null, null, "scan=1", 0, null, null, 0, null, DissociationType.Unknown, 1, null);

            PeptideSpectralMatch[] globalPsms = new PeptideSpectralMatch[1];
            Ms2ScanWithSpecificMass[] arrayOfSortedMS2Scans = { new Ms2ScanWithSpecificMass(scan, 0, 0, null) };
            CommonParameters CommonParameters = new CommonParameters(
                productMassTolerance: new PpmTolerance(5),
                scoreCutoff: 1,
                digestionParams: new DigestionParams(
                    maxMissedCleavages: 0,
                    minPeptideLength: 1,
                    maxModificationIsoforms: int.MaxValue,
                    initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain));
            ClassicSearchEngine cse = new ClassicSearchEngine(globalPsms, arrayOfSortedMS2Scans, new List<ModificationWithMass>(), new List<ModificationWithMass>(), new List<Protein> { prot }, new List<ProductType> { ProductType.B, ProductType.Y }, new OpenSearchMode(), CommonParameters, new List<string>());

            cse.Run();
            Assert.Less(globalPsms[0].Score, 2);
            Assert.Greater(globalPsms[0].Score, 1);
        }

        [Test]
        public static void TestVeryCloseExperimentalsModern()
        {
            IDictionary<int, List<Modification>> mods = new Dictionary<int, List<Modification>>();
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif);
            var prot = new Protein("MMMM", null, null, null, mods);
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 1);
            var thePep = prot.Digest(digestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass>()).First();

            var massArray = thePep.CompactPeptide(TerminusType.None).ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { ProductType.B, ProductType.Y });
            Array.Sort(massArray);
            double[] intensities = new double[] { 1, 1, 1, 1 };
            double[] mz = new double[] { 1, 2, massArray[4].ToMz(1), massArray[4].ToMz(1) + 1e-9 };
            MzSpectrum massSpectrum = new MzSpectrum(mz, intensities, false);
            MsDataScan scan = new MsDataScan(massSpectrum, 1, 1, true, Polarity.Positive, 1, new MzRange(300, 2000), "", MZAnalyzerType.Unknown, massSpectrum.SumOfAllY, null, null, "scan=1", 0, null, null, 0, null, DissociationType.Unknown, 1, null);

            PeptideSpectralMatch[] globalPsms = new PeptideSpectralMatch[1];
            Ms2ScanWithSpecificMass[] arrayOfSortedMS2Scans = { new Ms2ScanWithSpecificMass(scan, 600, 1, null) };
            CommonParameters CommonParameters = new CommonParameters(productMassTolerance: new PpmTolerance(5), scoreCutoff: 1, digestionParams: new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 1, maxModificationIsoforms: int.MaxValue, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain));

            var indexEngine = new IndexingEngine(new List<Protein> { prot }, new List<ModificationWithMass>(), new List<ModificationWithMass>(), new List<ProductType>
            { ProductType.B, ProductType.Y }, 1, DecoyType.Reverse, new List<DigestionParams> { CommonParameters.DigestionParams }, CommonParameters, 30000, new List<string>());
            var indexResults = (IndexingResults)indexEngine.Run();
            var cse = new ModernSearchEngine(globalPsms, arrayOfSortedMS2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, new List<ProductType> { ProductType.B, ProductType.Y }, 0, CommonParameters, new OpenSearchMode(), 0, new List<string>());

            cse.Run();
            Assert.Less(globalPsms[0].Score, 2);
            Assert.Greater(globalPsms[0].Score, 1);
        }

        [Test]
        public static void TestAllNaN()
        {
            IDictionary<int, List<Modification>> mods = new Dictionary<int, List<Modification>>();
            var prot = new Protein("XMMM", null, null, null, mods);
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 1);
            var thePep = prot.Digest(digestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass>()).First();

            var massArray = thePep.CompactPeptide(TerminusType.None).ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { ProductType.B });
            Array.Sort(massArray);
            double[] intensities = new double[] { 1, 1, 1, 1 };
            double[] mz = new double[] { 1, 2, 3, 4 };
            MzSpectrum massSpectrum = new MzSpectrum(mz, intensities, false);
            MsDataScan scan = new MsDataScan(massSpectrum, 1, 1, true, Polarity.Positive, 1, new MzRange(300, 2000), "", MZAnalyzerType.Unknown, massSpectrum.SumOfAllY, null, null, "scan=1", 0, null, null, 0, null, DissociationType.Unknown, 1, null);

            PeptideSpectralMatch[] globalPsms = new PeptideSpectralMatch[1];
            Ms2ScanWithSpecificMass[] arrayOfSortedMS2Scans = { new Ms2ScanWithSpecificMass(scan, 0, 0, null) };
            CommonParameters CommonParameters = new CommonParameters(productMassTolerance: new PpmTolerance(5), scoreCutoff: 1, digestionParams: new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 1, maxModificationIsoforms: int.MaxValue, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain));

            ClassicSearchEngine cse = new ClassicSearchEngine(globalPsms, arrayOfSortedMS2Scans, new List<ModificationWithMass>(), new List<ModificationWithMass>(), new List<Protein> { prot }, new List<ProductType> { ProductType.B, ProductType.Y }, new OpenSearchMode(), CommonParameters, new List<string>());

            cse.Run();
            Assert.IsNull(globalPsms[0]);
        }

        //[Test] MOVED TO MZLIB
        //public static void TestGoodPeptide()
        //public static void TestNoCleavage()
        //public static void TestBadPeptide()
        //public static void TestPeptideWithSetModifications()
        //public static void TestPeptideWithFixedModifications()
        //public static void TestDigestIndices()
        //public static void TestDigestDecoy()
        //public static void TestGoodPeptideWithLength()
    }
}
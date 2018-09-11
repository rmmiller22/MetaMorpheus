﻿using Chemistry;
using EngineLayer;
using EngineLayer.CrosslinkSearch;
using EngineLayer.Indexing;
using MassSpectrometry;
using Nett;
using NUnit.Framework;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.AminoAcidPolymer;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;
using UsefulProteomicsDatabases;
using OxyPlot;
using OxyPlot.Axes;
using OxyPlot.Series;
using OxyPlot.Annotations;
using MzLibUtil;

namespace Test
{
    [TestFixture]
    public static class XLTest
    {
        static IndexingResults indexResults { get; set; }
        static CommonParameters commonParameters { get; set; }
        static XlSearchParameters xlSearchParameters { get; set; }
        static List<Protein> proteinList { get; set; }
        static List<ModificationWithMass> variableModifications { get; set; }
        static List<ModificationWithMass> fixedModifications { get; set; }
        static List<ProductType> lp { get; set; }
        static CrosslinkerTypeClass crosslinker { get; set; }
        static List<PeptideWithSetModifications> digestedList { get; set; }

        private static void Generate()
        {
            //Generate parameters
            commonParameters = new CommonParameters(doPrecursorDeconvolution: false, cIons: true, zDotIons: true, scoreCutoff: 2, digestionParams: new DigestionParams(minPeptideLength: 5));
            xlSearchParameters = new XlSearchParameters { };

            //Create databases contain two protein.
            proteinList = new List<Protein> { new Protein("EKVLTSSAR", "Fake01"), new Protein("LSQKFPK", "Fake02") };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            ModificationWithMass mod1 = new ModificationWithMass("Oxidation of M", "Common Variable", motif1, TerminusLocalization.Any, 15.99491461957);
            ModificationMotif.TryGetMotif("C", out ModificationMotif motif2);
            ModificationWithMass mod2 = new ModificationWithMass("Carbamidomethyl of C", "Common Fixed", motif2, TerminusLocalization.Any, 57.02146372068994);
            variableModifications = new List<ModificationWithMass>() { mod1 };
            fixedModifications = new List<ModificationWithMass>() { mod2 };
            var localizeableModifications = new List<ModificationWithMass>();

            lp = new List<ProductType> { ProductType.BnoB1ions, ProductType.Y, ProductType.C, ProductType.Zdot };
            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();
            foreach (var mod in fixedModifications)
                modsDictionary.Add(mod, 0);
            int i = 1;
            foreach (var mod in variableModifications)
            {
                modsDictionary.Add(mod, (ushort)i);
                i++;
            }
            foreach (var mod in localizeableModifications)
            {
                modsDictionary.Add(mod, (ushort)i);
                i++;
            }

            //Generate digested peptide lists.
            digestedList = new List<PeptideWithSetModifications>();
            foreach (var item in proteinList)
            {
                var digested = item.Digest(commonParameters.DigestionParams, fixedModifications, variableModifications).ToList();
                digestedList.AddRange(digested);
            }
            foreach (var fdfd in digestedList)
            {
                fdfd.CompactPeptide(TerminusType.None);
            }

            //Run index engine
            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, lp, 1, DecoyType.Reverse, new List<DigestionParams>
            { commonParameters.DigestionParams }, commonParameters, 30000, new List<string>());
            indexResults = (IndexingResults)indexEngine.Run();

            //Generate crosslinker, which is DSSO here.
            crosslinker = new CrosslinkerTypeClass();
            crosslinker.SelectCrosslinker(xlSearchParameters.CrosslinkerType);
        }

        [Test]
        public static void XlTest_XlPosCal()
        {
            var prot = new Protein("MNNNKQQQQ", null);
            Protease protease = new Protease("New Custom Protease", new List<Tuple<string, FragmentationTerminus>> { new Tuple<string, FragmentationTerminus>("K", FragmentationTerminus.C) }, new List<Tuple<string, FragmentationTerminus>>(), CleavageSpecificity.Full, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            DigestionParams digestionParams = new DigestionParams(protease: protease.Name, minPeptideLength: 1, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            List<Modification> variableModifications = new List<Modification>();

            var ye = prot.Digest(digestionParams, new List<Modification>(), variableModifications).ToList();

            var pep = ye[0];
            Assert.AreEqual(pep.BaseSequence, "MNNNK");
            Crosslinker crosslinker = new Crosslinker();
            crosslinker.SelectCrosslinker(CrosslinkerType.DSS);
            Assert.AreEqual(crosslinker.CrosslinkerModSites, "K");
            Assert.AreEqual(Residue.GetResidue(crosslinker.CrosslinkerModSites).MonoisotopicMass, 128.09496301518999, 1e-9);
            var n = pep.Fragment(DissociationType.HCD, FragmentationTerminus.N);
            var c = pep.Fragment(DissociationType.HCD, FragmentationTerminus.C);
            Assert.AreEqual(n.Count(), 4);
            Assert.AreEqual(c.Count(), 4);
            Assert.AreEqual(c.First().NeutralMass, 146.10552769899999, 1e-6);
            var x = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(crosslinker.CrosslinkerModSites, pep).ToArray();
            Assert.AreEqual(x[0], 5);

            var pep2 = ye[2];
            Assert.AreEqual("MNNNKQQQQ", pep2.BaseSequence);
            var n2 = pep2.Fragment(DissociationType.HCD, FragmentationTerminus.N);
            var c2 = pep2.Fragment(DissociationType.HCD, FragmentationTerminus.C);
            Assert.AreEqual(n2.Count(), 8);
            Assert.AreEqual(c2.Count(), 8);
            var x2 = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(crosslinker.CrosslinkerModSites, pep2).ToArray();
            Assert.AreEqual(x2[0], 5);

            //Test crosslinker with multiple types of mod
            var protSTC = new Protein("GASTACK", null);
            var peps = protSTC.Digest(digestionParams, new List<Modification>(), variableModifications).ToList();
            var pepSTC = peps[0];
            Assert.AreEqual(pepSTC.BaseSequence, "GASTACK");
            Crosslinker crosslinker2 = new Crosslinker("ST", "C", "crosslinkerSTC", false, -18.01056, 0, 0, 0, 0, 0, 0);
            string crosslinkerModSitesAll = new string((crosslinker2.CrosslinkerModSites + crosslinker2.CrosslinkerModSites2).ToCharArray().Distinct().ToArray());
            Assert.AreEqual(crosslinkerModSitesAll, "STC");
        }

        [Test]
        public static void DeadendPeptideTest()
        {
            double[] intensity = new double[] { 1.1, 1.1, 0.5, 3.2, 0.5, 6.0 };
            int[] rank = CrosslinkSpectralMatch.GenerateIntensityRanks(intensity);
            int[] Rank = new int[] { 4, 3, 6, 2, 5, 1 };
            Assert.AreEqual(rank, Rank);
        }

        [Test]
        public static void XlTest_BSA_DSSO()
        {
            //Generate parameters
            var commonParameters = new CommonParameters(doPrecursorDeconvolution: false, dissociationType: DissociationType.EThcD, scoreCutoff: 1, digestionParams: new DigestionParams(minPeptideLength: 5));

            var xlSearchParameters = new XlSearchParameters { XlCharge_2_3_PrimeFragment = true };

            //Create databases contain two protein.
            var proteinList = new List<Protein> { new Protein("EKVLTSSAR", "Fake01"), new Protein("LSQKFPK", "Fake02") };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            Modification mod1 = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif1, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
            ModificationMotif.TryGetMotif("C", out ModificationMotif motif2);
            Modification mod2 = new Modification(_originalId: "Carbamidomethyl of C", _modificationType: "Common Fixed", _target: motif2, _locationRestriction: "Anywhere.", _monoisotopicMass: 57.02146372068994);
            var variableModifications = new List<Modification>() { mod1 };
            var fixedModifications = new List<Modification>() { mod2 };
            var localizeableModifications = new List<Modification>();

            var lp = new List<ProductType> { ProductType.b, ProductType.y, ProductType.c, ProductType.zPlusOne };
            Dictionary<Modification, ushort> modsDictionary = new Dictionary<Modification, ushort>();
            foreach (var mod in fixedModifications)
                modsDictionary.Add(mod, 0);
            int i = 1;
            foreach (var mod in variableModifications)
            {
                XlSearchParameters = new XlSearchParameters()
                {
                    XlPrecusorMsTl = new PpmTolerance(112000),
                    XlQuench_Tris = false,
                    XlQuench_H2O = false,
                    XlQuench_NH2 = true
                }
            };

            xLSearchTask.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabaseXl, false) }, new List<string> { myFileXl }, "test");
            xLSearchTask2.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabaseXl, false) }, new List<string> { myFileXl }, "test");
        }

            //Run index engine
            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType.Reverse, new List<DigestionParams>
            { commonParameters.DigestionParams }, commonParameters, 30000, new List<string>());

            var indexResults = (IndexingResults)indexEngine.Run();

            var fragmentIndexCount = indexResults.FragmentIndex.Count(p => p != null);
            var fragmentIndexAll = indexResults.FragmentIndex.Select((s, j) => new { j, s }).Where(p => p.s != null).Select(t => t.j).ToList();
            Assert.IsTrue(fragmentIndexAll.Count() > 0);

            //Get MS2 scans.
            var myMsDataFile = new XLTestDataFile();
            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, commonParameters.DoPrecursorDeconvolution, commonParameters.UseProvidedPrecursorInfo, commonParameters.DeconvolutionIntensityRatio, commonParameters.DeconvolutionMaxAssumedChargeState, commonParameters.DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

            //Generate crosslinker, which is DSSO here.
            Crosslinker crosslinker = new Crosslinker();
            crosslinker.SelectCrosslinker(xlSearchParameters.CrosslinkerType);

            //TwoPassCrosslinkSearchEngine.Run().
            List<CrosslinkSpectralMatch> newPsms = new List<CrosslinkSpectralMatch>();
            new TwoPassCrosslinkSearchEngine(newPsms, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, lp, 0, commonParameters, false, xlSearchParameters.XlPrecusorMsTl, crosslinker, xlSearchParameters.CrosslinkSearchTop, xlSearchParameters.CrosslinkSearchTopNum, xlSearchParameters.XlQuench_H2O, xlSearchParameters.XlQuench_NH2, xlSearchParameters.XlQuench_Tris, xlSearchParameters.XlCharge_2_3, xlSearchParameters.XlCharge_2_3_PrimeFragment, new List<string> { }).Run();
            
            foreach (var item in newPsms)
            {
                item.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0, 0, false);
            }

            //Test newPsms
            Assert.AreEqual(newPsms.Count(), 3);
            //Test Output
            var task = new XLSearchTask();
            task.WritePepXML_xl(newPsms, proteinList, null, variableModifications, fixedModifications, null, TestContext.CurrentContext.TestDirectory, "pep.XML", new List<string> { });
            task.WriteSingleToTsv(newPsms.Where(p => p.CrossType == PsmCrossType.Single).ToList(), TestContext.CurrentContext.TestDirectory, "singlePsms", new List<string> { });

            //Test PsmCross.XlCalculateTotalProductMasses
            var psmCrossAlpha = new CrosslinkSpectralMatch(digestedList[1], 0, 0, i, listOfSortedms2Scans[0], commonParameters.DigestionParams, new List<MatchedFragmentIon>());
            var psmCrossBeta = new CrosslinkSpectralMatch(digestedList[2], 0, 0, i, listOfSortedms2Scans[0], commonParameters.DigestionParams, new List<MatchedFragmentIon>());
            var linkPos = CrosslinkSpectralMatch.GetPossibleCrosslinkerModSites(crosslinker.CrosslinkerModSites, digestedList[1]);
            var productMassesAlphaList = CrosslinkedPeptide.XlGetTheoreticalFragmentIons(DissociationType.EThcD, false, crosslinker, linkPos, digestedList[2].MonoisotopicMass, digestedList[1]);
            Assert.AreEqual(productMassesAlphaList.First().Value.Count, 50); //TO DO: The number here should be manually verified.
        }

        [Test]
        public static void XlTest_BSA_DSS_file()
        {
            var task = Toml.ReadFile<XLSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData/XLSearchTaskconfig_BSA_DSS_23747.toml"), MetaMorpheusTask.tomlConfig);

            DbForTask db = new DbForTask(Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData/BSA.fasta"), false);
            string raw = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData/BSA_DSS_23747.mzML");
            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", task) }, new List<string> { raw }, new List<DbForTask> { db }, Path.Combine(Environment.CurrentDirectory, @"XlTestData")).Run();
        }

        [Test]
        public static void XlTest_GenerateUserDefinedCrosslinker()
        {
            XlSearchParameters xlSearchParameters = new XlSearchParameters();
            xlSearchParameters.CrosslinkerType = CrosslinkerType.UserDefined;
            xlSearchParameters.UdXLkerName = "CrossST-C";
            xlSearchParameters.UdXLkerResidues = "ST";
            xlSearchParameters.UdXLkerResidues2 = "C";
            xlSearchParameters.UdXLkerTotalMass = -18.01056;
            var crosslinker = XLSearchTask.GenerateUserDefinedCrosslinker(xlSearchParameters);
        }

        [Test]
        public static void XlTest_DiffCrosslinkSites()
        {
            //Generate parameters
            var commonParameters = new CommonParameters(doPrecursorDeconvolution: false, scoreCutoff: 1, digestionParams: new DigestionParams(minPeptideLength: 4));

            var xlSearchParameters = new XlSearchParameters
            {
                CrosslinkerType = CrosslinkerType.UserDefined,
                UdXLkerName = "CrossST-C",
                UdXLkerResidues = "ST",
                UdXLkerResidues2 = "C",
                UdXLkerTotalMass = -18.01056
            };

            //Create databases contain two protein.
            var proteinList = new List<Protein> { new Protein("VLTAR", "Fake01"), new Protein("LCQK", "Fake02") };

            ModificationMotif.TryGetMotif("M", out ModificationMotif motif1);
            Modification mod1 = new Modification(_originalId: "Oxidation of M", _modificationType: "Common Variable", _target: motif1, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491461957);
            var variableModifications = new List<Modification>() { mod1 };
            var fixedModifications = new List<Modification>();
            var localizeableModifications = new List<Modification>();

            var lp = new List<ProductType> { ProductType.b, ProductType.y };
            Dictionary<Modification, ushort> modsDictionary = new Dictionary<Modification, ushort>();

            int i = 1;
            foreach (var mod in variableModifications)
            {
                modsDictionary.Add(mod, (ushort)i);
                i++;
            }
            foreach (var mod in localizeableModifications)
            {
                modsDictionary.Add(mod, (ushort)i);
                i++;
            }

            //Generate digested peptide lists.
            List<PeptideWithSetModifications> digestedList = new List<PeptideWithSetModifications>();
            foreach (var item in proteinList)
            {
                var digested = item.Digest(commonParameters.DigestionParams, fixedModifications, variableModifications).ToList();
                digestedList.AddRange(digested);
            }

            //Run index engine
            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType.Reverse, new List<DigestionParams> { commonParameters.DigestionParams }, commonParameters, 30000, new List<string>());

            var indexResults = (IndexingResults)indexEngine.Run();

            //Get MS2 scans.
            var myMsDataFile = new XLTestDataFileDiffSite();
            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, commonParameters.DoPrecursorDeconvolution, commonParameters.UseProvidedPrecursorInfo, commonParameters.DeconvolutionIntensityRatio, commonParameters.DeconvolutionMaxAssumedChargeState, commonParameters.DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

            //Generate crosslinker, which is UserDefined here.
            var crosslinker = XLSearchTask.GenerateUserDefinedCrosslinker(xlSearchParameters);

            //TwoPassCrosslinkSearchEngine.Run().
            List<CrosslinkSpectralMatch> newPsms = new List<CrosslinkSpectralMatch>();
            new TwoPassCrosslinkSearchEngine(newPsms, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, lp, 0, commonParameters, false, xlSearchParameters.XlPrecusorMsTl, crosslinker, xlSearchParameters.CrosslinkSearchTop, xlSearchParameters.CrosslinkSearchTopNum, xlSearchParameters.XlQuench_H2O, xlSearchParameters.XlQuench_NH2, xlSearchParameters.XlQuench_Tris, xlSearchParameters.XlCharge_2_3, xlSearchParameters.XlCharge_2_3_PrimeFragment, new List<string> { }).Run();
            Assert.AreEqual(newPsms.Count(), 1);
        }

        /// <summary>
        /// Verifies that crosslinker is generated properly
        /// </summary>
        [Test]
        public static void CrosslinkCreateTest()
        {
            Assert.That((XLSearchTask.GenerateUserDefinedCrosslinker(new XlSearchParameters())).GetType().Equals(typeof(CrosslinkerTypeClass)));
        }

        /// <summary>
        /// Makes sure helper methods that generate indices function properly
        /// </summary>
        [Test]
        public static void XLSearchWithGeneratedIndices()
        {
            XLSearchTask xlSearchTask = new XLSearchTask();
            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData\BSA_DSSO_ETchD6010.mgf");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData\BSA.fasta");
            string folderPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestXLSearch");
            DbForTask db = new DbForTask(myDatabase, false);
            List<(string, MetaMorpheusTask)> taskList = new List<(string, MetaMorpheusTask)> { ("TestXLSearch", xlSearchTask) };

            //creates .params files if they do not exist
            xlSearchTask.RunTask(Path.Combine(folderPath, @"CreateParams"), new List<DbForTask> { db }, new List<string> { myFile }, "normal");
            //tests .params files
            xlSearchTask.RunTask(Path.Combine(folderPath, @"TestParams"), new List<DbForTask> { db }, new List<string> { myFile }, "normal");

            var baseDir = Path.GetDirectoryName(db.FilePath);
            var directory = new DirectoryInfo(baseDir);
            DirectoryInfo[] directories = directory.GetDirectories();
            foreach (DirectoryInfo possibleFolder in directories)
            {
                if (File.Exists(Path.Combine(possibleFolder.FullName, "indexEngine.params")))
                {
                    File.Delete(possibleFolder.GetFiles().ElementAt(0).FullName);
                }
            }
            //tests without .params files
            xlSearchTask.RunTask(Path.Combine(folderPath, @"TestNoParams"), new List<DbForTask> { db }, new List<string> { myFile }, "normal");

            var lines = File.ReadAllLines(Path.Combine(folderPath, @"CreateParams\xl_intra_fdr.tsv"));
            var lines2 = File.ReadAllLines(Path.Combine(folderPath, @"TestParams\xl_intra_fdr.tsv"));
            var lines3 = File.ReadAllLines(Path.Combine(folderPath, @"TestNoParams\xl_intra_fdr.tsv"));

            Assert.That(lines.SequenceEqual(lines2) && lines2.SequenceEqual(lines3));
            
        }

        [Test]
        public static void Xl_TestDataFilePSMControl()
        {
            Generate();

            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData/BSA_DSSO_ETchD6010.mgf");
            CommonParameters commonParameters = new CommonParameters(doPrecursorDeconvolution: false);
            MyFileManager myFileManager = new MyFileManager(true);
            var msDataFile = myFileManager.LoadFile(filePath, 300, 0.01, true, true, commonParameters);
            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(msDataFile, filePath, commonParameters.DoPrecursorDeconvolution, commonParameters.UseProvidedPrecursorInfo, commonParameters.DeconvolutionIntensityRatio, commonParameters.DeconvolutionMaxAssumedChargeState, commonParameters.DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

            //TwoPassCrosslinkSearchEngine.Run().
            List<PsmCross> newPsms = new List<PsmCross>();
            new TwoPassCrosslinkSearchEngine(newPsms, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, lp, 0, commonParameters, false, false, false, xlSearchParameters.XlPrecusorMsTl, crosslinker, xlSearchParameters.CrosslinkSearchTop, xlSearchParameters.CrosslinkSearchTopNum, xlSearchParameters.XlQuench_H2O, xlSearchParameters.XlQuench_NH2, xlSearchParameters.XlQuench_Tris, xlSearchParameters.XlCharge_2_3, new List<string> { }).Run();

            var compactPeptideToProteinPeptideMatch = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>();

            new CrosslinkAnalysisEngine(newPsms, compactPeptideToProteinPeptideMatch, proteinList, variableModifications, fixedModifications, lp, null, crosslinker, TerminusType.None, commonParameters, new List<string> { }).Run();
            foreach (var item in newPsms)
            {
                item.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0, 0, false);
            }

            DrawPeptideSpectralMatch(msDataFile.GetAllScansList()[0], newPsms.First());
        }

        private static Dictionary<ProductType, OxyColor> productTypeDrawColors = new Dictionary<ProductType, OxyColor>
        { { ProductType.B, OxyColors.Blue },
          { ProductType.BnoB1ions, OxyColors.Blue },
          { ProductType.Y, OxyColors.Red },
          { ProductType.C, OxyColors.ForestGreen },
          { ProductType.Zdot, OxyColors.Orange },
         { ProductType.X, OxyColors.Violet },
        { ProductType.None, OxyColors.Purple }};

        public static void DrawPeptideSpectralMatch(MsDataScan msDataScan, PsmCross psmToDraw)
        {
            // x is m/z, y is intensity
            var spectrumMzs = msDataScan.MassSpectrum.XArray;
            var spectrumIntensities = msDataScan.MassSpectrum.YArray;

            string subTitle = psmToDraw.FullSequence;
            if (psmToDraw.BetaPsmCross != null)
            {
                subTitle += "-" + psmToDraw.BetaPsmCross.FullSequence;
            }
            PlotModel model = new PlotModel { Title = "Spectrum Annotation of Scan #" + msDataScan.OneBasedScanNumber, DefaultFontSize = 15, Subtitle = subTitle };
            model.Axes.Add(new LinearAxis { Position = AxisPosition.Bottom, Title = "m/z", Minimum = 0, Maximum = spectrumMzs.Max() * 1.2, AbsoluteMinimum = 0, AbsoluteMaximum = spectrumMzs.Max() * 5 });
            model.Axes.Add(new LinearAxis { Position = AxisPosition.Left, Title = "Intensity", Minimum = 0, Maximum = spectrumIntensities.Max() * 1.2, AbsoluteMinimum = 0, AbsoluteMaximum = spectrumIntensities.Max() * 1.3 });
            model.Axes[1].Zoom(0, spectrumIntensities.Max() * 1.1);

            LineSeries[] allIons = new LineSeries[spectrumMzs.Length];

            // draw the matched peaks; if the PSM is null, we're just drawing the peaks in the scan without annotation, so skip this part
            if (psmToDraw != null)
            {
                if (psmToDraw.BetaPsmCross != null)
                {
                    psmToDraw.MatchedIons = psmToDraw.MatchedIons.Concat(psmToDraw.BetaPsmCross.MatchedIons).ToList();
                }
                foreach (var peak in psmToDraw.MatchedIons)
                {
                    OxyColor ionColor = productTypeDrawColors[peak.TheoreticalFragmentIon.ProductType];

                    int i = msDataScan.MassSpectrum.GetClosestPeakIndex(peak.Mz).Value;

                    // peak line
                    allIons[i] = new LineSeries();
                    allIons[i].Color = ionColor;
                    allIons[i].StrokeThickness = 1;
                    allIons[i].Points.Add(new DataPoint(peak.Mz, 0));
                    allIons[i].Points.Add(new DataPoint(peak.Mz, spectrumIntensities[i]));

                    // peak annotation
                    var peakAnnotation = new TextAnnotation();
                    peakAnnotation.TextRotation = -60;
                    peakAnnotation.Font = "Arial";
                    peakAnnotation.FontSize = 4;
                    peakAnnotation.FontWeight = 1.0;
                    peakAnnotation.TextColor = ionColor;
                    peakAnnotation.StrokeThickness = 0;
                    //string gly = peak.TheoreticalFragmentIon.ProductType == ProductType.None ? "-Hex" : "";
                    peakAnnotation.Text = "(" + peak.Mz.ToString("F3") + "@+" + peak.TheoreticalFragmentIon.Charge.ToString() + ") " + peak.TheoreticalFragmentIon.ProductType.ToString().ToLower().First() + "-" + peak.TheoreticalFragmentIon.IonNumber;
                    peakAnnotation.TextPosition = new DataPoint(allIons[i].Points[1].X, allIons[i].Points[1].Y + peakAnnotation.Text.Length * 1.5 / 4);
                    peakAnnotation.TextHorizontalAlignment = HorizontalAlignment.Left;
                    model.Annotations.Add(peakAnnotation);

                    model.Series.Add(allIons[i]);
                }
            }

            for (int i = 0; i < spectrumMzs.Length; i++)
            {
                allIons[i] = new LineSeries();
                allIons[i].Color = OxyColors.DimGray;
                allIons[i].StrokeThickness = 0.5;
                allIons[i].Points.Add(new DataPoint(spectrumMzs[i], 0));
                allIons[i].Points.Add(new DataPoint(spectrumMzs[i], spectrumIntensities[i]));
                model.Series.Add(allIons[i]);
            }

            // Axes are created automatically if they are not defined

            // Set the Model property, the INotifyPropertyChanged event will make the WPF Plot control update its content
            using (var stream = File.Create(Path.Combine(TestContext.CurrentContext.TestDirectory, psmToDraw.ScanNumber + ".pdf")))
            {
                PdfExporter pdf = new PdfExporter { Width = 500, Height = 210 };
                pdf.Export(model, stream);
            }
        }
    }

    internal class XLTestDataFile : MsDataFile
    {
        //Create DSSO crosslinked fake MS data. Include single, deadend, loop, inter, intra crosslinks ms2 data for match.
        public XLTestDataFile() : base(2, new SourceFile(null, null, null, null, null))
        {
            var mz1 = new double[] { 1994.05.ToMz(3), 846.4963.ToMz(1), 1004.495.ToMz(1), 1093.544.ToMz(1), 1043.561.ToMz(1) };
            var intensities1 = new double[] { 1, 1, 1, 1, 1 };
            var MassSpectrum1 = new MzSpectrum(mz1, intensities1, false);
            var ScansHere = new List<MsDataScan> { new MsDataScan(MassSpectrum1, 1, 1, true, Polarity.Positive, 1, new MzLibUtil.MzRange(0, 10000), "ff", MZAnalyzerType.Unknown, 1000, 1, null, "scan=1") };

            var mz2 = new double[] { 100, 201.1234, 244.1656, 391.2340, 420.2201, 521.2678, 634.3519, 889.965, 1044.568, 1094.551, 1279.671, 1378.74, 1491.824 };
            var intensities2 = new double[] { 100, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
            var MassSpectrum2 = new MzSpectrum(mz2, intensities2, false);
            ScansHere.Add(new MsDataScan(MassSpectrum2, 2, 2, true, Polarity.Positive, 1.0,
                new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Unknown, 112, 1.0, null, "scan=2", 1994.05.ToMz(3),
                3, 1, 1994.05.ToMz(3), 2, DissociationType.HCD, 1, 1994.05.ToMz(3)));

            var mz3 = new double[] { 100, 201.1234, 244.1656, 391.2340 };
            var intensities3 = new double[] { 100, 1, 1, 1 };
            var MassSpectrum3 = new MzSpectrum(mz3, intensities3, false);
            ScansHere.Add(new MsDataScan(MassSpectrum3, 3, 2, true, Polarity.Positive, 1.0,
                new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Unknown, 103, 1.0, null, "scan=3", 846.4963.ToMz(1),
                1, 1, 846.4963.ToMz(1), 2, DissociationType.HCD, 1, 846.4963.ToMz(1)));

            var mz4 = new double[] { 100, 201.1234, 244.1656, 391.2340 };
            var intensities4 = new double[] { 100, 1, 1, 1 };
            var MassSpectrum4 = new MzSpectrum(mz4, intensities4, false);
            ScansHere.Add(new MsDataScan(MassSpectrum4, 4, 2, true, Polarity.Positive, 1.0,
                new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Unknown, 103, 1.0, null, "scan=4", 1004.491.ToMz(1),
                1, 1, 1004.491.ToMz(1), 2, DissociationType.HCD, 1, 1004.491.ToMz(1)));

            Scans = ScansHere.ToArray();
        }

        public string FilePath
        {
            get
            {
                return "XLTestDataFile";
            }
        }

        public string Name
        {
            get
            {
                return "XLTestDataFile";
            }
        }

        public void ReplaceFirstScanArrays(double[] mz, double[] intensities)
        {
            MzSpectrum massSpectrum = new MzSpectrum(mz, intensities, false);
            Scans[0] = new MsDataScan(massSpectrum, Scans[0].OneBasedScanNumber, Scans[0].MsnOrder, Scans[0].IsCentroid, Scans[0].Polarity, Scans[0].RetentionTime, Scans[0].ScanWindowRange, Scans[0].ScanFilter, Scans[0].MzAnalyzer, massSpectrum.SumOfAllY, Scans[0].InjectionTime, null, Scans[0].NativeId);
        }
    }

    internal class XLTestDataFileDiffSite : MsDataFile
    {
        //Create DSSO crosslinked fake MS data. Include single, deadend, loop, inter, intra crosslinks ms2 data for match.
        public XLTestDataFileDiffSite() : base(2, new SourceFile(null, null, null, null, null))
        {
            var mz1 = new double[] { 100, 1030.5956.ToMz(1) };
            var intensities1 = new double[] { 100, 1 };
            var MassSpectrum1 = new MzSpectrum(mz1, intensities1, false);
            var ScansHere = new List<MsDataScan> { new MsDataScan(MassSpectrum1, 1, 1, true, Polarity.Positive, 1, new MzLibUtil.MzRange(0, 10000), "ff", MZAnalyzerType.Unknown, 1000, 1, null, "scan=1") };

            var mz2 = new double[] { 100, 147.1128, 175.119, 213.1598, 246.1561, 275.1714, 757.4388, 786.4541, 819.4504, 857.4912, 885.4974, 918.5189, 932.5345 };
            var intensities2 = new double[] { 100, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
            var MassSpectrum2 = new MzSpectrum(mz2, intensities2, false);
            ScansHere.Add(new MsDataScan(MassSpectrum2, 2, 2, true, Polarity.Positive, 1.0,
                new MzLibUtil.MzRange(0, 10000), "f", MZAnalyzerType.Unknown, 112, 1.0, null, "scan=2", 1030.5956.ToMz(1),
                1, 1, 1030.5956.ToMz(1), 2, DissociationType.HCD, 1, 1030.5956.ToMz(1)));

            Scans = ScansHere.ToArray();
        }

        public string FilePath
        {
            get
            {
                return "XLTestDataFileDiffSite";
            }
        }

        public string Name
        {
            get
            {
                return "XLTestDataFileDiffSite";
            }
        }

        public void ReplaceFirstScanArrays(double[] mz, double[] intensities)
        {
            MzSpectrum massSpectrum = new MzSpectrum(mz, intensities, false);
            Scans[0] = new MsDataScan(massSpectrum, Scans[0].OneBasedScanNumber, Scans[0].MsnOrder, Scans[0].IsCentroid, Scans[0].Polarity, Scans[0].RetentionTime, Scans[0].ScanWindowRange, Scans[0].ScanFilter, Scans[0].MzAnalyzer, massSpectrum.SumOfAllY, Scans[0].InjectionTime, null, Scans[0].NativeId);
        }
    }
}
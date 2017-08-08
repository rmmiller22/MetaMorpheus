﻿using Proteomics;
using System;
using System.Linq;
using System.Collections.Generic;

namespace EngineLayer
{
    [Serializable]
    public class CompactPeptide : CompactPeptideBase
    {
        #region Public Constructors

        public CompactPeptide(PeptideWithSetModifications peptideWithSetModifications)
        {
            double theMass = 0;

            if (peptideWithSetModifications.allModsOneIsNterminus.TryGetValue(1, out ModificationWithMass pep_n_term_variable_mod))
                foreach (double nl in pep_n_term_variable_mod.neutralLosses)
                    theMass = pep_n_term_variable_mod.monoisotopicMass - nl;
            else
                theMass = 0;
            NTerminalMasses = ComputeFollowingFragmentMasses(peptideWithSetModifications, theMass, 1, 1).ToArray();

            theMass = 0;
            if (peptideWithSetModifications.allModsOneIsNterminus.TryGetValue(peptideWithSetModifications.Length + 2, out ModificationWithMass pep_c_term_variable_mod))
                foreach (double nl in pep_c_term_variable_mod.neutralLosses)
                    theMass = pep_c_term_variable_mod.monoisotopicMass - nl;
            else
                theMass = 0;
            CTerminalMasses = ComputeFollowingFragmentMasses(peptideWithSetModifications, theMass, peptideWithSetModifications.Length, -1).ToArray();

            MonoisotopicMassIncludingFixedMods = peptideWithSetModifications.MonoisotopicMass;
        }

        public CompactPeptide(PeptideWithSetModifications peptideWithSetModifications, TerminusType terminusType)
        {
            double theMass = 0;
            NTerminalMasses = new double[1];
            CTerminalMasses = new double[1];
            if (terminusType==TerminusType.None || terminusType == TerminusType.N)
            {
                if (peptideWithSetModifications.allModsOneIsNterminus.TryGetValue(1, out ModificationWithMass pep_n_term_variable_mod))
                    foreach (double nl in pep_n_term_variable_mod.neutralLosses)
                        theMass = pep_n_term_variable_mod.monoisotopicMass - nl;
                else
                    theMass = 0;
                NTerminalMasses = ComputeFollowingFragmentMasses(peptideWithSetModifications, theMass, 1, 1).ToArray();
            }
            if (terminusType == TerminusType.None || terminusType == TerminusType.C)
            {
                theMass = 0;
                if (peptideWithSetModifications.allModsOneIsNterminus.TryGetValue(peptideWithSetModifications.Length + 2, out ModificationWithMass pep_c_term_variable_mod))
                    foreach (double nl in pep_c_term_variable_mod.neutralLosses)
                        theMass = pep_c_term_variable_mod.monoisotopicMass - nl;
                else
                    theMass = 0;
                CTerminalMasses = ComputeFollowingFragmentMasses(peptideWithSetModifications, theMass, peptideWithSetModifications.Length, -1).ToArray();
            }
            MonoisotopicMassIncludingFixedMods = peptideWithSetModifications.MonoisotopicMass;
        }

        #endregion Public Constructors

        #region Public Methods

        public double[] ProductMassesMightHaveDuplicatesAndNaNs(List<ProductType> productTypes)
        {
            int massLen = 0;
            bool containsAdot = productTypes.Contains(ProductType.Adot);
            bool containsB = productTypes.Contains(ProductType.B);
            bool containsBnoB1 = productTypes.Contains(ProductType.BnoB1ions);
            bool containsC = productTypes.Contains(ProductType.C);
            bool containsX = productTypes.Contains(ProductType.X);
            bool containsY = productTypes.Contains(ProductType.Y);
            bool containsZdot = productTypes.Contains(ProductType.Zdot);

            if (containsAdot)
                throw new NotImplementedException();
            if (containsBnoB1)
                massLen += NTerminalMasses.Length - 1;
            if (containsB)
                massLen += NTerminalMasses.Length;
            if (containsC)
                massLen += NTerminalMasses.Length;
            if (containsX)
                throw new NotImplementedException();
            if (containsY)
                massLen += CTerminalMasses.Length;
            if (containsZdot)
                massLen += CTerminalMasses.Length;

            double[] massesToReturn = new double[massLen];

            int i = 0;
            for (int j = 0; j < NTerminalMasses.Length; j++)
            {
                var hm = NTerminalMasses[j];
                if (containsBnoB1 && j > 0)
                {
                    massesToReturn[i] = hm;
                    i++;
                }
                if (containsB)
                {
                    massesToReturn[i] = hm;
                    i++;
                }
                if (containsC)
                {
                    massesToReturn[i] = hm + nitrogenAtomMonoisotopicMass + 3 * hydrogenAtomMonoisotopicMass;
                    i++;
                }
            }
            for (int j = 0; j < CTerminalMasses.Length; j++)
            {
                var hm = CTerminalMasses[j];
                if (containsY)
                {
                    massesToReturn[i] = hm + waterMonoisotopicMass;
                    i++;
                }
                if (containsZdot)
                {
                    massesToReturn[i] = hm + oxygenAtomMonoisotopicMass - nitrogenAtomMonoisotopicMass;
                    i++;
                }
            }
            return massesToReturn;
        }

        #endregion Public Methods

    }
}
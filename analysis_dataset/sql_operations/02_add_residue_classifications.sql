-- ============================================================================
-- Migration: 002_add_residue_classifications.sql
-- Description: Add residue type classification columns
-- 
-- This migration adds two types of classifications:
-- 1. Chemical property classification (restype): hydrophobic, aromatic, polar, etc.
-- 2. Molecule type classification (type): protein, chlorophyll, carotenoid
-- 
-- These classifications enable analysis by residue chemistry and molecule class
-- ============================================================================

-- ----------------------------------------------------------------------------
-- SECTION 1: Chemical Property Classification (restype_i)
-- Purpose: Classify first residue by amino acid chemical properties
-- ----------------------------------------------------------------------------

ALTER TABLE data ADD COLUMN restype_i TEXT;

UPDATE data
SET restype_i = CASE resname_i
    -- Hydrophobic (aliphatic) amino acids
    WHEN 'ALA' THEN 'hydrophobic'
    WHEN 'VAL' THEN 'hydrophobic'
    WHEN 'LEU' THEN 'hydrophobic'
    WHEN 'ILE' THEN 'hydrophobic'
    WHEN 'MET' THEN 'hydrophobic'
    
    -- Aromatic amino acids
    WHEN 'PHE' THEN 'aromatic'
    WHEN 'TYR' THEN 'aromatic'
    WHEN 'TRP' THEN 'aromatic'
    
    -- Polar uncharged amino acids
    WHEN 'SER' THEN 'polar'
    WHEN 'THR' THEN 'polar'
    WHEN 'ASN' THEN 'polar'
    WHEN 'GLN' THEN 'polar'
    
    -- Negatively charged (acidic) amino acids
    WHEN 'ASP' THEN 'ASP'
    WHEN 'GLU' THEN 'GLU'
    
    -- Positively charged (basic) amino acids
    WHEN 'LYS' THEN 'basic'
    WHEN 'ARG' THEN 'basic'
    WHEN 'HIS' THEN 'basic'
    
    -- Special cases (unique structural properties)
    WHEN 'GLY' THEN 'special'  -- No side chain
    WHEN 'PRO' THEN 'special'  -- Cyclic structure
    WHEN 'CYS' THEN 'special'  -- Disulfide bonds
    
    -- Chlorophyll cofactors
    WHEN 'CLA' THEN 'chlorophyll'
    WHEN 'CLB' THEN 'chlorophyll'
    WHEN 'CHL' THEN 'chlorophyll'
    
    -- Carotenoid cofactors
    WHEN 'LUT' THEN 'carotenoid'  -- Lutein
    WHEN 'VIO' THEN 'carotenoid'  -- Violaxanthin
    WHEN 'XAT' THEN 'carotenoid'  -- Xanthophyll
    WHEN 'NEO' THEN 'carotenoid'  -- Neoxanthin
    WHEN 'BCR' THEN 'carotenoid'  -- Beta-carotene
    WHEN 'NEX' THEN 'carotenoid'  -- Neoxanthin
    
    ELSE 'unknown'
END;


-- ----------------------------------------------------------------------------
-- SECTION 2: Chemical Property Classification (restype_j)
-- Purpose: Classify second residue by amino acid chemical properties
-- ----------------------------------------------------------------------------

ALTER TABLE data ADD COLUMN restype_j TEXT;

UPDATE data
SET restype_j = CASE resname_j
    -- Hydrophobic (aliphatic) amino acids
    WHEN 'ALA' THEN 'hydrophobic'
    WHEN 'VAL' THEN 'hydrophobic'
    WHEN 'LEU' THEN 'hydrophobic'
    WHEN 'ILE' THEN 'hydrophobic'
    WHEN 'MET' THEN 'hydrophobic'
    
    -- Aromatic amino acids
    WHEN 'PHE' THEN 'aromatic'
    WHEN 'TYR' THEN 'aromatic'
    WHEN 'TRP' THEN 'aromatic'
    
    -- Polar uncharged amino acids
    WHEN 'SER' THEN 'polar'
    WHEN 'THR' THEN 'polar'
    WHEN 'ASN' THEN 'polar'
    WHEN 'GLN' THEN 'polar'
    
    -- Negatively charged (acidic) amino acids
    WHEN 'ASP' THEN 'ASP'
    WHEN 'GLU' THEN 'GLU'
    
    -- Positively charged (basic) amino acids
    WHEN 'LYS' THEN 'basic'
    WHEN 'ARG' THEN 'basic'
    WHEN 'HIS' THEN 'basic'
    
    -- Special cases (unique structural properties)
    WHEN 'GLY' THEN 'special'  -- No side chain
    WHEN 'PRO' THEN 'special'  -- Cyclic structure
    WHEN 'CYS' THEN 'special'  -- Disulfide bonds
    
    -- Chlorophyll cofactors
    WHEN 'CLA' THEN 'chlorophyll'
    WHEN 'CLB' THEN 'chlorophyll'
    WHEN 'CHL' THEN 'chlorophyll'
    
    -- Carotenoid cofactors
    WHEN 'LUT' THEN 'carotenoid'  -- Lutein
    WHEN 'VIO' THEN 'carotenoid'  -- Violaxanthin
    WHEN 'XAT' THEN 'carotenoid'  -- Xanthophyll
    WHEN 'NEO' THEN 'carotenoid'  -- Neoxanthin
    WHEN 'BCR' THEN 'carotenoid'  -- Beta-carotene
    WHEN 'NEX' THEN 'carotenoid'  -- Neoxanthin
    
    ELSE 'unknown'
END;


-- ----------------------------------------------------------------------------
-- SECTION 3: Molecule Type Classification (type_i)
-- Purpose: Classify first residue by molecule type (protein vs cofactor)
-- ----------------------------------------------------------------------------

ALTER TABLE data ADD COLUMN type_i TEXT;

UPDATE data
SET type_i = CASE resname_i
    -- All 20 standard amino acids are classified as 'protein'
    WHEN 'ALA' THEN 'protein'
    WHEN 'VAL' THEN 'protein'
    WHEN 'LEU' THEN 'protein'
    WHEN 'ILE' THEN 'protein'
    WHEN 'MET' THEN 'protein'
    WHEN 'PHE' THEN 'protein'
    WHEN 'TYR' THEN 'protein'
    WHEN 'TRP' THEN 'protein'
    WHEN 'SER' THEN 'protein'
    WHEN 'THR' THEN 'protein'
    WHEN 'ASN' THEN 'protein'
    WHEN 'GLN' THEN 'protein'
    WHEN 'ASP' THEN 'protein'
    WHEN 'GLU' THEN 'protein'
    WHEN 'LYS' THEN 'protein'
    WHEN 'ARG' THEN 'protein'
    WHEN 'HIS' THEN 'protein'
    WHEN 'GLY' THEN 'protein'
    WHEN 'PRO' THEN 'protein'
    WHEN 'CYS' THEN 'protein'
    
    -- Chlorophyll cofactors (grouped as 'chlorophyll')
    WHEN 'CLA' THEN 'chlorophyll'
    WHEN 'CLB' THEN 'chlorophyll'
    WHEN 'CHL' THEN 'chlorophyll'
    
    -- Carotenoid cofactors (kept as specific types for analysis)
    WHEN 'LUT' THEN 'LUT'  -- Lutein
    WHEN 'VIO' THEN 'VIO'  -- Violaxanthin
    WHEN 'XAT' THEN 'XAT'  -- Xanthophyll
    WHEN 'NEO' THEN 'NEO'  -- Neoxanthin
    WHEN 'BCR' THEN 'BCR'  -- Beta-carotene
    WHEN 'NEX' THEN 'NEO'  -- Neoxanthin (same as NEO)
    
    ELSE 'unknown'
END;


-- ----------------------------------------------------------------------------
-- SECTION 4: Molecule Type Classification (type_j)
-- Purpose: Classify second residue by molecule type (protein vs cofactor)
-- ----------------------------------------------------------------------------

ALTER TABLE data ADD COLUMN type_j TEXT;

UPDATE data
SET type_j = CASE resname_j
    -- All 20 standard amino acids are classified as 'protein'
    WHEN 'ALA' THEN 'protein'
    WHEN 'VAL' THEN 'protein'
    WHEN 'LEU' THEN 'protein'
    WHEN 'ILE' THEN 'protein'
    WHEN 'MET' THEN 'protein'
    WHEN 'PHE' THEN 'protein'
    WHEN 'TYR' THEN 'protein'
    WHEN 'TRP' THEN 'protein'
    WHEN 'SER' THEN 'protein'
    WHEN 'THR' THEN 'protein'
    WHEN 'ASN' THEN 'protein'
    WHEN 'GLN' THEN 'protein'
    WHEN 'ASP' THEN 'protein'
    WHEN 'GLU' THEN 'protein'
    WHEN 'LYS' THEN 'protein'
    WHEN 'ARG' THEN 'protein'
    WHEN 'HIS' THEN 'protein'
    WHEN 'GLY' THEN 'protein'
    WHEN 'PRO' THEN 'protein'
    WHEN 'CYS' THEN 'protein'
    
    -- Chlorophyll cofactors (grouped as 'chlorophyll')
    WHEN 'CLA' THEN 'chlorophyll'
    WHEN 'CLB' THEN 'chlorophyll'
    WHEN 'CHL' THEN 'chlorophyll'
    
    -- Carotenoid cofactors (kept as specific types for analysis)
    WHEN 'LUT' THEN 'LUT'  -- Lutein
    WHEN 'VIO' THEN 'VIO'  -- Violaxanthin
    WHEN 'XAT' THEN 'XAT'  -- Xanthophyll
    WHEN 'NEO' THEN 'NEO'  -- Neoxanthin
    WHEN 'BCR' THEN 'BCR'  -- Beta-carotene
    WHEN 'NEX' THEN 'NEO'  -- Neoxanthin (same as NEO)
    
    ELSE 'unknown'
END;

-- ============================================================================
-- End of migration 002
-- ============================================================================


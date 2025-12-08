-- ============================================================================
-- Migration: 001_add_psbs_modifiers.sql
-- Description: Add derived columns for data transformation and identification
-- 
-- This migration adds:
-- 1. Concatenated residue identifiers for easier reference
-- 2. Unit conversions for lifetime data
-- 3. PsbS-specific residue renumbering and chain labeling
-- ============================================================================

-- ----------------------------------------------------------------------------
-- SECTION 1: Concatenated Residue Identifiers
-- Purpose: Create human-readable identifiers combining residue ID and name
-- ----------------------------------------------------------------------------

-- Add resid_resname_i column (e.g., "123ALA")
-- Combines residue ID and name for the first residue in the pair
ALTER TABLE data ADD COLUMN resid_resname_i TEXT;

UPDATE data
SET resid_resname_i = CAST(resid_i AS TEXT) || resname_i;

-- Add resid_resname_j column (e.g., "456GLY")
-- Combines residue ID and name for the second residue in the pair
ALTER TABLE data ADD COLUMN resid_resname_j TEXT;

UPDATE data
SET resid_resname_j = CAST(resid_j AS TEXT) || resname_j;

-- Add resid_resname_pair column (e.g., "123ALA-456GLY")
-- Combines both residue identifiers for the complete interaction pair
ALTER TABLE data ADD COLUMN resid_resname_pair TEXT;

UPDATE data
SET resid_resname_pair = resid_resname_i || '-' || resid_resname_j;


-- ----------------------------------------------------------------------------
-- SECTION 2: Unit Conversions
-- Purpose: Convert lifetime from nanoseconds to microseconds for convenience
-- ----------------------------------------------------------------------------

-- Add lifetime_microseconds column
-- Converts lifetime_ns to microseconds (1 microsecond = 1000 nanoseconds)
ALTER TABLE data ADD COLUMN lifetime_microseconds REAL;

UPDATE data
SET lifetime_microseconds = lifetime_ns / 1000.0;


-- ----------------------------------------------------------------------------
-- SECTION 3: PsbS-Specific Modifiers
-- Purpose: Renumber PsbS residues and create chain labels
-- Context: PsbS has 212 residues offset that needs correction
-- ----------------------------------------------------------------------------

-- Add resid_j_renumbered column
-- Renumbers PsbS residue IDs by subtracting 212 for values greater than 212
-- This corrects for sequence offset in the PsbS structure
ALTER TABLE data ADD COLUMN resid_j_renumbered INTEGER;

UPDATE data
SET resid_j_renumbered = CASE 
    WHEN resid_j > 212 THEN resid_j - 212
    ELSE resid_j
END;

-- Add chain_j_renamed column
-- Assigns chain labels based on residue position:
-- - 'A' for residues 1-212 (first domain)
-- - 'B' for residues >212 (second domain)
ALTER TABLE data ADD COLUMN chain_j_renamed TEXT;

UPDATE data
SET chain_j_renamed = CASE 
    WHEN resid_j > 212 THEN 'B'
    ELSE 'A'
END;

-- Add resid_j_renumbered_resname column (e.g., "23GLY" for renumbered residue)
-- Combines the renumbered residue ID with its name
ALTER TABLE data ADD COLUMN resid_j_renumbered_resname TEXT;

UPDATE data
SET resid_j_renumbered_resname = CAST(resid_j_renumbered AS TEXT) || resname_j;

-- ============================================================================
-- End of migration 001
-- ============================================================================



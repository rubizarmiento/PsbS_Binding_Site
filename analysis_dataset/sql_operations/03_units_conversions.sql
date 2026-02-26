-- ============================================================================
-- Migration: 001_add_psbs_modifiers.sql
-- Description: Add derived columns for data transformation and identification
-- 
-- This migration adds:
-- 1. Unit conversions for lifetime data
-- ============================================================================
-- ----------------------------------------------------------------------------
-- SECTION 1: Unit Conversions
-- Purpose: Convert lifetime from nanoseconds to microseconds for convenience
-- ----------------------------------------------------------------------------

-- Adds the columns sum_ms, median_ms, and p90_ms to store the converted values in milliseconds
-- Converts lifetime_ns to microseconds (1 microsecond = 1000 nanoseconds)
ALTER TABLE data ADD COLUMN sum_ms REAL;
ALTER TABLE data ADD COLUMN median_ms REAL;
ALTER TABLE data ADD COLUMN p90_ms REAL;


UPDATE data
SET sum_ms = sum_ns / 1000.0,
    median_ms = median_ns / 1000.0,
    p90_ms = p90_ns / 1000.0;   




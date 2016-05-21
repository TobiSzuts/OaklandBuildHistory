# -*- coding: utf-8 -*-
"""
Created on Sat May 21 13:55:34 2016

@author: tobis
"""

# OaklandBuildHistory

This project explores the year houses were built in Oakland, focused on several neighborhoods (not the entire city, at least to begin with).

Initially, it will display on a map the year each house was built.  This will be a chloropleth map, with the polygons being the parcel numbers.

Basic strategy is to use parcel information from the county of Alameda to find each parcel in a neighborhood (by distance from a point), then use reverse geocode API from Google to find the street number, then query property info for that address on Zillow.


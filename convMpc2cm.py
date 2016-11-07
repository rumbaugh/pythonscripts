import math as m
import sys

def convMpc2cm(d='1',inverse=False,abb=False):
    #converts Mpc to cm. If inverse is true, it's the other way around.
    #if abb=True, it gives the answer in units of 10^24 cm
    if ((inverse) & (abb)): sys.exit("Abbreviation not compatible with inverse conversion")
    conversion = 3.085677581e24
    if inverse: conversion = 1.0/conversion
    if abb: conversion = 3.085677581
    newd = d*conversion
    return newd

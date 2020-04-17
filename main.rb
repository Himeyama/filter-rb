#!/usr/bin/env ruby
require "csv"

require "./test"

data = []
CSV.foreach("test.csv").with_index do |row, i|
    data[i] = row.map(&:to_f)
end

data = NArray[*data]
t = data[0, true]
x = data[1, true]

s = 1925
fp = 40
fs = 59
gpass = 1
gstop = 40
fn = s / 2.0
wp = fp / fn
ws = fs / fn

signal = Filter.new

n, wn = signal.buttord(wp, ws, gpass, gstop)

p n, wn
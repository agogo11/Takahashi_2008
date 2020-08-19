;\\ Find the array index for which the array value is closest to the desired value.
FUNCTION yzclosest, in_array, in_value, $
  noabsolute = noabsolute, $
  difference=difference, $
  value=value

  diff = abs(in_array - in_value)
  IF KEYWORD_SET(noabsolute) THEN diff = in_array - in_value
  pt = (WHERE(diff EQ MIN(diff)))[0]
  difference = diff[pt]
  value = in_array[pt]
  RETURN, pt

END
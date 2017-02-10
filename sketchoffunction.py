def convertTGASfile(filename):
  '''this function will:
      -read in a TGAS FITS file
      -convert it into a table
      -extract the RA, Dec, parallax, etc... from the table
      -trim to only the good data
      -calculate distance, x, y, z for every star
      -create some demonstration plots (labeled by the filename)
        +x vs y
        +y vs z
        +ra vs dec
      -save these plots to files (whose names depend on the input "filename")
      -save a text file that contains X, Y, Z positions and aboslute magnitudes
      

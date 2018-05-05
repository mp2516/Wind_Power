"""1) IMPORTING LIBRARIES AND READING NETCDF FILES"""

#importing the appropriate libraries
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset as NetCDFfile
 
# reading the netCDF files ('r' is for read only)
# reading each file into a different name, this program reads the files for 1 year only
nc_1979_01 = NetCDFfile('C:\\Users\\18072\\OneDrive\\Documents\\Physics\\1st_Year\\Summer Project\\Data\\wnd10m.gdas.197901.grb2.nc', mode ='r')
nc_1979_02 = NetCDFfile('C:\\Users\\18072\\OneDrive\\Documents\\Physics\\1st_Year\\Summer Project\\Data\\wnd10m.gdas.197902.grb2.nc', mode ='r')
nc_1979_03 = NetCDFfile('C:\\Users\\18072\\OneDrive\\Documents\\Physics\\1st_Year\\Summer Project\\Data\\wnd10m.gdas.197903.grb2.nc', mode ='r')
nc_1979_04 = NetCDFfile('C:\\Users\\18072\\OneDrive\\Documents\\Physics\\1st_Year\\Summer Project\\Data\\wnd10m.gdas.197904.grb2.nc', mode ='r')
nc_1979_05 = NetCDFfile('C:\\Users\\18072\\OneDrive\\Documents\\Physics\\1st_Year\\Summer Project\\Data\\wnd10m.gdas.197905.grb2.nc', mode ='r')
nc_1979_06 = NetCDFfile('C:\\Users\\18072\\OneDrive\\Documents\\Physics\\1st_Year\\Summer Project\\Data\\wnd10m.gdas.197906.grb2.nc', mode ='r')
nc_1979_07 = NetCDFfile('C:\\Users\\18072\\OneDrive\\Documents\\Physics\\1st_Year\\Summer Project\\Data\\wnd10m.gdas.197907.grb2.nc', mode ='r')
nc_1979_08 = NetCDFfile('C:\\Users\\18072\\OneDrive\\Documents\\Physics\\1st_Year\\Summer Project\\Data\\wnd10m.gdas.197908.grb2.nc', mode ='r')
nc_1979_09 = NetCDFfile('C:\\Users\\18072\\OneDrive\\Documents\\Physics\\1st_Year\\Summer Project\\Data\\wnd10m.gdas.197909.grb2.nc', mode ='r')
nc_1979_10 = NetCDFfile('C:\\Users\\18072\\OneDrive\\Documents\\Physics\\1st_Year\\Summer Project\\Data\\wnd10m.gdas.197910.grb2.nc', mode ='r')
nc_1979_11 = NetCDFfile('C:\\Users\\18072\\OneDrive\\Documents\\Physics\\1st_Year\\Summer Project\\Data\\wnd10m.gdas.197911.grb2.nc', mode ='r')
nc_1979_12 = NetCDFfile('C:\\Users\\18072\\OneDrive\\Documents\\Physics\\1st_Year\\Summer Project\\Data\\wnd10m.gdas.197912.grb2.nc', mode ='r')

# opening the land file
# the land file is 0 everywhere where there is sea and 1 on land
# this file removes all the points
nc_land = NetCDFfile('C:\\Users\\18072\\OneDrive\\Documents\\Physics\\1st Year\\Summer Project\\flxf01.gdas.1990.grb2.nc', mode ='r')
land = np.array(nc_land.variables['LAND_L1_Avg'][1])
nc_land.close()


u01 = np.array(nc_1979_01.variables['U_GRD_L103'][:])
v01 = np.array(nc_1979_01.variables['V_GRD_L103'][:])
u02 = np.array(nc_1979_02.variables['U_GRD_L103'][:])
v02 = np.array(nc_1979_02.variables['V_GRD_L103'][:])
u03 = np.array(nc_1979_03.variables['U_GRD_L103'][:])
v03 = np.array(nc_1979_03.variables['V_GRD_L103'][:])
u04 = np.array(nc_1979_04.variables['U_GRD_L103'][:])
v04 = np.array(nc_1979_04.variables['V_GRD_L103'][:])
u05 = np.array(nc_1979_05.variables['U_GRD_L103'][:])
v05 = np.array(nc_1979_05.variables['V_GRD_L103'][:])
u06 = np.array(nc_1979_06.variables['U_GRD_L103'][:])
v06 = np.array(nc_1979_06.variables['V_GRD_L103'][:])
u07 = np.array(nc_1979_07.variables['U_GRD_L103'][:])
v07 = np.array(nc_1979_07.variables['V_GRD_L103'][:])
u08 = np.array(nc_1979_08.variables['U_GRD_L103'][:])
v08 = np.array(nc_1979_08.variables['V_GRD_L103'][:])
u09 = np.array(nc_1979_09.variables['U_GRD_L103'][:])
v09 = np.array(nc_1979_09.variables['V_GRD_L103'][:])
u10 = np.array(nc_1979_10.variables['U_GRD_L103'][:])
v10 = np.array(nc_1979_10.variables['V_GRD_L103'][:])
u11 = np.array(nc_1979_11.variables['U_GRD_L103'][:])
v11 = np.array(nc_1979_11.variables['V_GRD_L103'][:])
u12 = np.array(nc_1979_12.variables['U_GRD_L103'][:])
v12 = np.array(nc_1979_12.variables['V_GRD_L103'][:])

# now defining arrays
# u = , v = , lat = latitude, lon = longitude, sp = speed
# u and v arrays are concatenated along the time axis
# u and v have dimensions of (35,48,8760) as they are 365 x 24 = 8760 in the time axis.
u = np.concatenate((u01, u02, u03, u04, u05, u06, u07, u08, u09, u10, u11, u12), axis=0)
v = np.concatenate((v01, v02, v03, v04, v05, v06, v07, v08, v09, v10, v11, v12), axis=0)
lat = np.array(nc_1979_01.variables['lat'][:])
lon = np.array(nc_1979_01.variables['lon'][:])
sp = np.sqrt(u**2 + v**2)

# close the file now that I have finished
nc_1979_01.close()
nc_1979_02.close()
nc_1979_03.close()
nc_1979_04.close()
nc_1979_05.close()
nc_1979_06.close()
nc_1979_07.close()
nc_1979_08.close()
nc_1979_09.close()
nc_1979_10.close()
nc_1979_11.close()
nc_1979_12.close()

###############################################

"""2) DEFINING FUNCTIONS"""

def size_square(lat,lon):
    
    #defining the Earth's radius
    R_earth = 6.371*(10**6)
    
    """Longitude"""
    #shape = 48, start: -11.8750, end: 2.8125
    #average difference = 0.3060
    #difference between consecutive values in an array
    lon_diff = np.ediff1d(lon)
    #in this case the difference is constant so this line is redundant
    mean_lon_diff = np.mean(lon_diff)
    #length of arc of a circle is radius*angle
    #where angle is in radians
    lon_dis = (R_earth)*((mean_lon_diff*(np.pi))/180)
    
    """Latitude"""
    #shape = 35, start: 59.7918, end: 49.1760
    #average difference = 0.3033
    #same calculations as before
    #reversing array before starting to ensure the difference is postive
    lat_diff = np.ediff1d(lat[::-1])
    mean_lat_diff = np.mean(lat_diff)
    #the only difference is the multiplication by cos(lat)
    #this is because the distance representated by a constant
    #latitude change becomes smaller as latitude increases
    #remember that python works in radians not degrees for cos
    lat_dis = (np.cos((lat*(np.pi))/180))*(R_earth)*((mean_lat_diff*(np.pi))/180)
    
    """Area"""
    #shape = 35
    # the area of the square is now the multiplication of distances
    # reshaped so as to obey broadcasting rules further down
    square = (lat_dis * lon_dis).reshape(1,35,1)
    return square
    
    
# Two general power arrays #

def power():
    
    """Power per Turbine"""
    # define the regions for the function
    cutinspeed = sp<=3.5
    #the use of logical_and is because the 'truth value of an array is ambiguous'
    vcubed = np.logical_and(3.5 <= sp, sp < 13)
    maxpower = np.logical_and(13 <= sp, sp < 25)
    cutoffspeed = sp>=25
    
    # numerically evaluate the values of each point
    # P is zero when x<3.5m/s or when x>25m/s
    P = np.where(cutinspeed, 0.0, 0.0)
    # P has a cubic relationship when 3.5 < x < 13m/s
    # coefficients are calculated using boundary condtions
    P = np.where(vcubed, 1393*(sp**3) - 59711, P)
    # P has a constant max power of 3MW (Vestas V90 3MW turbine) when 13 < x < 25m/s
    P = np.where(maxpower, 3*(10**6), P)
    P = np.where(cutoffspeed, 0.0, P)
    power_turbine = P
    
    """Power per Square"""
    # space per wind turbine is 1km squared
    spacing = 1.0e6
    
    # number of turbines in one square = square area / area per turbine
    turbines_square = (size_square(lat,lon))/(spacing)
    #dimensions of P_turbine are (720, 35, 48)
    #so Turbines_square must be reshaped to (1,1,48) to obey broadcasting rules
    power_per_square = (power_turbine)*(turbines_square)
    #power_per_square = power_per_square[:,:,:16]
    
    """Land coverage"""
    # remove any points that are not on land as it was found that otherwise all the points appear to be in the Atlantic!
    power_per_square = power_per_square * land
    
    power_per_square = power_per_square[:,:,:]
    return power_per_square

  
def power_square(x,y):
    # power per square given x and y
    # Changing the arguments of the function so they are easier to use
    # Whether the time axis is sliced or not did not seem to make any difference
    # to the time it takes the code to run
    power_for_square = power()[:,x,y]
    return power_for_square


# Two mean power arrays #

def mean_power():
    # time averaged power without x and y
    # mean_power is (1,35,48)
    mean_power = np.mean(power(), axis=0)
    return mean_power
    
    
    
def mean_power_square(x,y):
    # time-averaged power for a square given x and y
    mean_power_square = mean_power()[x,y]
    return mean_power_square
    

# Max power #
    
def max_power():
    # finding the value of the maximum power across all points
    # this is the highest value in the average power array
    max_power = mean_power().max()
    return max_power
    


###############################################

"""3) CORRELATION"""

def correlation(P_time,x1,y1):
    # correlation coefficients between a point x1,y1 and a point x2,y2
    correlation = np.corrcoef(P_time, power_square(x1,y1))
    # the correlation coefficient produces a 2x2 matrix for every pair of points
    return correlation[0,1]


  
def correlation_coefficients(P_time):
    # Defining the empty arrays
    # Change the (35,48) to e.g. (10,11) if you want the code to run faster
    # (5,6) takes approximately 3 seconds
    cor_coef = np.zeros([35,48])
    
    # the correlation coefficients are calculated relative to the max power point

    # Change the range for faster code (see above)
    for x in range(35):
        for y in range(48):
            # producing the array of correlation coefficient for point (xpos,ypos)
            cor_coef[x,y] = correlation(P_time,x,y)
            # save the file correlation coefficients to the disc
            # the file takes a long time to create so this is for quicker testing
            cor_coef = np.nan_to_num(cor_coef)
            #np.save('correlation_coefficients', cor_coef)
    # the function if called will run the code again and overwrite any files
    return cor_coef

    

def best_points(num):
    # finding the index of the point with maximum power
    # defining an empty array
    P_total = 0.0
    P_time = 0.0
    power1 = power()[:,:,:]
    mean_power = np.mean(power1, axis=0)
    for i in range(num):
        xpoint, ypoint = np.where(mean_power == np.amax(mean_power))
        max_power_ix = int(xpoint[0])
        max_power_iy = int(ypoint[0])
        mean_power[max_power_ix, max_power_iy] = 0.0
        power1[:,max_power_ix, max_power_iy] = 0.0
        P_time += power_square(max_power_ix, max_power_iy)
        P_total += mean_power.max()
        max_power = mean_power.max()
        print('latitude = %-10g longitude = %-10g power = %-14g cumulative = %-10g' % \
          (lat[max_power_ix], lon[max_power_iy], max_power, P_total))

    # load the correlation coefficent file from the disc to save time
    #cor_coef = np.load('correlation_coefficients.npy')
    
    # P_point starts of at max_power and is built up until the demand of the country is met
    
    
    while P_total <= 4.1096e10 * 0.5:
        metric = (((correlation_coefficients(P_time)) * -1) + 1) * mean_power
        for i in range(5):
            # where function returns a tuple so requires two variables to be input into
            x_cor, y_cor = np.where((metric) == np.amax(metric))
            # find array elements only where land coverage = 1
            # converts the tuple into integers so we can index the longitude and latitude arrays
            xpoint = int(x_cor[0])
            ypoint = int(y_cor[0])
            # make the correlation coefficient for that point 0 so it is undesirable
            mean_power[(xpoint), (ypoint)] = 0
            power1[:,xpoint, ypoint] = 0.0
            metric[(xpoint), (ypoint)] = 0        
            # P_total is the running total of power
            P_total += mean_power_square(xpoint, ypoint)
            P_time += power_square(xpoint, ypoint)
            # the correlation coefficient varies between -1 and 1, when multiplied by -1 everything that was negatively /
            # correlated becomes positive so we can maximise the function
        
            # print in columns so that the information is readable
            print (' latitude = %-10g longitude = %-10g power = %-14g cumulative = %-10g' % \
                   (lat[xpoint], lon[ypoint], mean_power_square(xpoint, ypoint), P_total))
            # when the demand of the country is reached, exit the function
    
    t = np.arange(0, 8760, 1.0)
    power_quota = np.zeros(len(t))
    for i in range(8760):
        power_quota[i] = (4.1096e9)*(np.cos(((2*(np.pi))/8760)*(t[i]))) + 4.1096e10
   
    shortfall = power_quota - P_time
    count = 0
    for i in range(8760):
        if shortfall[i] < 0:
            count += 1
    print (count)
            
    
    while shortfall.any() > 0 and P_total < 4.1096e10:
        # Which is the most underperforming site?
        t1 = (np.where((shortfall) == np.amax(shortfall)))
        t_underperforming = int(t1[0])
        # Have co-ordinates of worst performing site.
            # Now zero the shortfall so this point isn't chosen again
        x, y = np.where(power1[t_underperforming,:,:] == np.amax(power1[t_underperforming,:,:]))
        xtime = int(x[0])
        ytime = int(y[0])
        power1[:,xtime, ytime] = 0.0
        P_time += power_square(xtime, ytime)
        P_total += mean_power_square(xtime, ytime)
        shortfall = power_quota - P_time
        print ('latitude = %-10g longitude = %-10g power = %-14g cumulative = %-10g' % \
               (lat[xtime], lon[ytime], mean_power_square(xtime, ytime), P_total))                           
    
    
    #np.save('P_time_150', P_time)
    return P_total, P_time


def plotting():
    t = np.arange(0, 8760, 1.0)
    #P_time = np.load('P_time.npy')
    P_total, P_time = best_points(30)
    y1 = P_time
    y2 = np.zeros(8760)
    for i in range(8760):
        y2[i] = (4.1096e9)*(np.cos(((2*(np.pi))/8760)*(t[i]))) + 4.1096e10

    shortfall = y2 - P_time
    count = 0
    storage = 0
    for i in range(8760):
        if shortfall[i] < 0:
            count += 1
            storage += shortfall[i]
    print (count)
    print (storage)
    
    plt.plot(t, y1, 'r', t, y2, 'b')
    plt.xlabel('Time in hours')
    plt.ylabel('Power')
    plt.legend(['Power provided by wind turbines', 'Power demand of the country'])
    plt.title('Power vs Time')
    plt.savefig('Power.pdf')
    plt.show()
    
plotting()

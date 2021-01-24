import numpy as np
import matplotlib.pyplot as plt

def change_degree_range(value):
    """ -> from -180 to 180"""
    value = np.where(value > 180, value - 360, value)
    value = np.where(value < -180, value + 360, value)
    return value

def view_mesh(x, y, filename):
    plt.plot(x, y, c="k")
    plt.plot(np.transpose(x), np.transpose(y), c="k")
    plt.axis('equal')
    #plt.show()
    plt.savefig(filename)
    plt.close()
    return


def merc(lat, lon):
    """メルカトル投影"""
    # degree -> rad
    lat = np.deg2rad(lat)
    lon = np.deg2rad(lon)

    x = lon
    y = np.arctanh(np.sin(lat))

    # rad -> deg
    x = np.rad2deg(x)
    y = np.rad2deg(y)
    return x, y

def spherical2cartesian(r, theta_rad, phi_rad):
    x = r*np.sin(theta_rad)*np.cos(phi_rad)
    y = r*np.sin(theta_rad)*np.sin(phi_rad)
    z = r*np.cos(theta_rad)
    return x, y, z

def cartesian2spherical(x, y, z):
    r = (x**2 + y**2 + z**2)**0.5
    theta_rad = np.arccos(z / (x**2 + y**2 + z**2)**0.5)
    phi_rad = np.where(y>=0, 1, -1)*np.arccos(x/(x**2 + y**2)**0.5)
    return r, theta_rad, phi_rad

def calc_new_pole(lat_point_1, lon_point_1, lat_point_2, lon_point_2):
    """
        Point1とPoint2を通る大円を新たな赤道とするような新たな極の座標を求める
    """
    # deg -> rad
    lat_point_1 = np.deg2rad(lat_point_1)
    lon_point_1 = np.deg2rad(lon_point_1)
    lat_point_2 = np.deg2rad(lat_point_2)
    lon_point_2 = np.deg2rad(lon_point_2)

    # 緯度を余緯度に変換（球面座標系の演算規則を使いたいので）
    alpha_1 = np.deg2rad(90) - lat_point_1
    alpha_2 = np.deg2rad(90) - lat_point_2

    # Point1とPoint2の外積をとる
    # (=原点とPoint1とPoint2を含む面に直交するベクトルが得られる)

    x_1, y_1, z_1 = spherical2cartesian(1, alpha_1, lon_point_1)
    x_2, y_2, z_2 = spherical2cartesian(1, alpha_2, lon_point_2)

    #外積
    x_pole = y_1*z_2 - z_1*y_2
    y_pole = z_1*x_2 - x_1*z_2
    z_pole = x_1*y_2 - y_1*x_2

    _, alpha_pole, lon_pole = cartesian2spherical(x_pole, y_pole, z_pole)

    # 余緯度→緯度へ変換
    lat_pole = np.deg2rad(90) - alpha_pole

    # rad -> deg
    lat_pole = np.rad2deg(lat_pole)
    lon_pole = np.rad2deg(lon_pole)
    return lat_pole, lon_pole

def oblique(lat, lon, lat_point_1, lon_point_1, lat_point_2, lon_point_2):
    """射軸変換
        Point1とPoint2を通る大円を新たな赤道とするような射軸変換を行う
        その後、中心の経度がcenter_lonになるように回転する
        lat, lon : 変換対象のLat, Lon
    """

    # 変換後、北極になる座標を求める
    pole_lat, pole_lon = calc_new_pole(lat_point_1, lon_point_1, lat_point_2, lon_point_2)

    # degree -> rad
    lat = np.deg2rad(lat)
    lon = np.deg2rad(lon)
    pole_lat = np.deg2rad(pole_lat)
    pole_lon = np.deg2rad(pole_lon)

    # 射軸変換
    # Lat
    v = np.sin(pole_lat)*np.sin(lat) + np.cos(pole_lat)*np.cos(lat)*np.cos(lon - pole_lon)
    v = np.where(np.isclose(v, 1), 1, v)
    v = np.where(np.isclose(v,-1), -1, v)
    lat_converted = np.arcsin(v)

    # Lon
    v = (np.sin(pole_lat)*np.sin(lat_converted) - np.sin(lat))/(np.cos(pole_lat)*np.cos(lat_converted))
    v = np.where(np.isclose(v, 1), 1, v)
    v = np.where(np.isclose(v,-1), -1, v)
    lon_converted = np.where(np.sin(lon - pole_lon)>=0, 1, -1)*np.fabs(np.arccos(v))

    # rad -> deg
    lat_converted = np.rad2deg(lat_converted)
    lon_converted = np.rad2deg(lon_converted) + 180 #中心経度を0度に表示
    return lat_converted, lon_converted

def main():

    ###########################################################################
    # SETTINGS
    ###########################################################################
    # mesh grid（描画範囲）
    lat, lon = np.meshgrid(np.linspace(10, 60, 11), np.linspace(100, 160, 11))

    # 射軸変換用の大円(=赤道)になる2点
    lat_point_1 = 34.0
    lon_point_1 = 140.0
    lat_point_2 = 34.0
    lon_point_2 = 120.0
    # この場合、中心軸のLonは130.0degree
    ###########################################################################
    # END OF SETTINGS
    ###########################################################################


    # 射軸変換
    lat, lon = oblique(lat, lon, lat_point_1, lon_point_1, lat_point_2, lon_point_2)

    # 定義域の調整
    lat = change_degree_range(lat)
    lon = change_degree_range(lon)

    # メルカトル投影
    x, y = merc(lat, lon)

    # 描画
    view_mesh(x, y, "oblique_merc.png")
    return

if __name__ == "__main__":
    main()
    pass

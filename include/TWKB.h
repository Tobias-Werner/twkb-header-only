//
// Created by Tobias Werner
//

#ifndef TWKB_HEADER_ONLY_TWKB_H
#define TWKB_HEADER_ONLY_TWKB_H

#include <vector>

namespace TWKB {

    using namespace std;

    using bytes_t = vector<unsigned char>;

    struct PosXY {
        double x, y;

        PosXY(double x, double y) : x(x), y(y) {}
    };

    struct PosXYZ {
        double x, y, z;

        PosXYZ(double x, double y, double z) : x(x), y(y), z(z) {}
    };

    struct PosXYZT {
        double x, y, z, t;

        PosXYZT(double x, double y, double z, double t) : x(x), y(y), z(z), t(t) {}
    };


    class GeomFactory {

    private:

        void setTypeAndPrecision(bytes_t &twkb, unsigned char type, char precisionXY) {
            auto &header = twkb[0];
            header |= type;
            header |= encodeZigZag(precisionXY) << 4;
        }

        void setMetadataBits(bytes_t &twkb, bool bbox, bool sizes, bool idlist, bool extendedDimensions, bool emptyGeom) {
            auto &metadataHeader = twkb[1];
            metadataHeader |= (bbox) ? 0x01 : 0x00;
            metadataHeader |= (sizes) ? 0x02 : 0x00;
            metadataHeader |= (idlist) ? 0x04 : 0x00;
            metadataHeader |= (extendedDimensions) ? 0x08 : 0x00;
            metadataHeader |= (emptyGeom) ? 0x10 : 0x00;
        }

        void addZDimensions(bytes_t &twkb, unsigned char precisionZ) {
            unsigned char extendedInformation = 0x01;
            extendedInformation |= (precisionZ << 2);
            twkb.push_back(extendedInformation);
        }

        void addZTDimensions(bytes_t &twkb, unsigned char precisionZ, unsigned char precisionT) {
            unsigned char extendedInformation = 0x03;
            extendedInformation |= (precisionZ << 2);
            extendedInformation |= (precisionT << 5);
            twkb.push_back(extendedInformation);
        }

        int shrink(double value, char precision) {
            return (precision >= 0) ? round(value * pow10(precision)) : round(value / pow10(-precision));
        }

        bytes_t encode(double value, char precision) {
            return encodeVarint(encodeZigZag(shrink(value, precision)));
        }

        void append(bytes_t &twkb) {}

        template<typename T, typename... Types>
        void append(bytes_t &twkb, T var1, Types... var2) {
            twkb.insert(twkb.end(), var1.begin(), var1.end());
            append(twkb, var2...);
        }

    public:

        static int pow10(unsigned char base) {
            switch (base) {
                case 0 :
                    return 1;
                case 1 :
                    return 10;
                case 2 :
                    return 100;
                case 3 :
                    return 1000;
                case 4 :
                    return 10000;
                case 5 :
                    return 100000;
                case 6 :
                    return 1000000;
                case 7 :
                    return 10000000;
                case 8 :
                    return 100000000;
                default:
                    return 0;
            }
        }

        static unsigned int decodeVarint(bytes_t bytes) {

            unsigned int result = 0x00;

            for (size_t i = bytes.size(); i > 0; i--) {
                result |= (bytes.at(i - 1) & 0x7F);
                result <<= (i - 1 > 0) ? 7 : 0;
            }

            return result;
        }

        static bytes_t encodeVarint(unsigned int value) {
            bytes_t result;

            do {
                result.push_back(0x80 | (value & 0x7F));
                value >>= 7;
            } while (value > 0x00);

            result.back() &= 0x7F;

            return result;
        }

        static unsigned int encodeZigZag(int value) {
            return (value << 1) ^ (value >> 31);
        }

        static int decodeZigZag(unsigned int value) {
            return -(value & 1) ^ (value >> 1);
        }

        bytes_t makePoint(PosXY location, char precisionXY) {

            bytes_t twkb({0x00, 0x00});

            // Set type and precision
            setTypeAndPrecision(twkb, 0x01, precisionXY);

            // Prepare metadata header
            // bbox | size | idlist | extended dimensions | empty geom
            setMetadataBits(twkb, false, false, false, false, false);

            // Locations
            auto x = encode(location.x, precisionXY);
            auto y = encode(location.y, precisionXY);

            append(twkb, x, y);

            return twkb;
        }

        bytes_t makePoint(PosXYZ location, char precisionXY, char precisionZ) {

            bytes_t twkb({0x00, 0x00});

            // Set type and precision
            setTypeAndPrecision(twkb, 0x01, precisionXY);

            // Prepare metadata header
            // bbox | size | idlist | extended dimensions | empty geom
            setMetadataBits(twkb, false, false, false, true, false);

            addZDimensions(twkb, precisionZ);

            auto x = encode(location.x, precisionXY);
            auto y = encode(location.y, precisionXY);
            auto z = encode(location.z, precisionZ);

            append(twkb, x, y, z);

            return twkb;
        }

        bytes_t makePoint(PosXYZT location, char precisionXY, char precisionZ, char precisionT) {

            bytes_t twkb({0x00, 0x00});

            // Set type and precision
            setTypeAndPrecision(twkb, 0x01, precisionXY);

            // Prepare metadata header
            // bbox | size | idlist | extended dimensions | empty geom
            setMetadataBits(twkb, false, false, false, true, false);

            addZTDimensions(twkb, precisionZ, precisionT);

            auto x = encode(location.x, precisionXY);
            auto y = encode(location.y, precisionXY);
            auto z = encode(location.z, precisionZ);
            auto t = encode(location.t, precisionT);

            append(twkb, x, y, z, t);

            return twkb;
        }

        bytes_t makeLine(vector<PosXY> locations, char precisionXY, bool bbox) {
            bytes_t twkb({0x00, 0x00});

            // Set type and precision
            setTypeAndPrecision(twkb, 0x02, precisionXY);

            // Prepare metadata header
            // bbox | size | idlist | extended dimensions | empty geom
            setMetadataBits(twkb, bbox, false, false, false, false);

            if (bbox) {

                int minx = shrink(locations.front().x, precisionXY);
                int maxx = minx;
                int miny = shrink(locations.front().y, precisionXY);
                int maxy = miny;

                for (auto &location : locations) {
                    int tmpX = shrink(location.x, precisionXY);
                    int tmpY = shrink(location.y, precisionXY);

                    if (tmpX < minx) minx = tmpX;
                    if (tmpY < miny) miny = tmpY;
                    if (tmpX > maxx) maxx = tmpX;
                    if (tmpY > maxy) maxy = tmpY;
                }

                bytes_t bboxX = encodeVarint(encodeZigZag(minx));
                bytes_t bboxY = encodeVarint(encodeZigZag(miny));
                bytes_t deltaX = encodeVarint(encodeZigZag(maxx - minx));
                bytes_t deltaY = encodeVarint(encodeZigZag(maxy - miny));

                append(twkb, bboxX, deltaX, bboxY, deltaY);

            }

            // Set number of containing points
            bytes_t npoints = encodeVarint(locations.size());

            append(twkb, npoints);

            int xShrinked = shrink(locations.front().x, precisionXY);
            int yShrinked = shrink(locations.front().y, precisionXY);

            bytes_t x = encodeVarint(encodeZigZag(xShrinked));
            bytes_t y = encodeVarint(encodeZigZag(yShrinked));

            append(twkb, x, y);

            for (int i = 1; i < locations.size(); i++) {

                int deltaX = shrink(locations[i].x, precisionXY) - xShrinked;
                int deltaY = shrink(locations[i].y, precisionXY) - yShrinked;

                x = encodeVarint(encodeZigZag(deltaX));
                y = encodeVarint(encodeZigZag(deltaY));

                append(twkb, x, y);

                xShrinked = shrink(locations[i].x, precisionXY);
                yShrinked = shrink(locations[i].y, precisionXY);

            }

            return twkb;
        }

        bytes_t makeLine(vector<PosXYZ> locations, char precisionXY, char precisionZ, bool bbox) {
            bytes_t twkb({0x00, 0x00});

            // Set type and precision
            setTypeAndPrecision(twkb, 0x02, precisionXY);

            // Prepare metadata header
            // bbox | size | idlist | extended dimensions | empty geom
            setMetadataBits(twkb, bbox, false, false, true, false);

            // Add extendet dimension information
            addZDimensions(twkb, precisionZ);

            // Add bounding box - if selected
            if (bbox) {

                int minx = shrink(locations.front().x, precisionXY);
                int maxx = minx;
                int miny = shrink(locations.front().y, precisionXY);
                int maxy = miny;
                int minz = shrink(locations.front().z, precisionZ);
                int maxz = minz;

                for (auto &location : locations) {
                    int tmpX = shrink(location.x, precisionXY);
                    int tmpY = shrink(location.y, precisionXY);
                    int tmpZ = shrink(location.z, precisionZ);

                    if (tmpX < minx) minx = tmpX;
                    if (tmpY < miny) miny = tmpY;
                    if (tmpZ < minz) minz = tmpZ;
                    if (tmpX > maxx) maxx = tmpX;
                    if (tmpY > maxy) maxy = tmpY;
                    if (tmpZ > maxz) maxz = tmpZ;
                }

                bytes_t bboxX = encodeVarint(encodeZigZag(minx));
                bytes_t bboxY = encodeVarint(encodeZigZag(miny));
                bytes_t bboxZ = encodeVarint(encodeZigZag(minz));
                bytes_t deltaX = encodeVarint(encodeZigZag(maxx - minx));
                bytes_t deltaY = encodeVarint(encodeZigZag(maxy - miny));
                bytes_t deltaZ = encodeVarint(encodeZigZag(maxz - minz));

                append(twkb, bboxX, deltaX, bboxY, deltaY, bboxZ, deltaZ);

            }


            // Set number of containing points
            bytes_t npoints = encodeVarint(locations.size());
            twkb.insert(twkb.end(), npoints.begin(), npoints.end());

            int xShrinked = shrink(locations.front().x, precisionXY);
            int yShrinked = shrink(locations.front().y, precisionXY);
            int zShrinked = shrink(locations.front().z, precisionZ);

            bytes_t x = encodeVarint(encodeZigZag(xShrinked));
            bytes_t y = encodeVarint(encodeZigZag(yShrinked));
            bytes_t z = encodeVarint(encodeZigZag(zShrinked));

            append(twkb, x, y, z);

            for (int i = 1; i < locations.size(); i++) {

                int deltaX = shrink(locations[i].x, precisionXY) - xShrinked;
                int deltaY = shrink(locations[i].y, precisionXY) - yShrinked;
                int deltaZ = shrink(locations[i].z, precisionZ) - zShrinked;

                x = encodeVarint(encodeZigZag(deltaX));
                y = encodeVarint(encodeZigZag(deltaY));
                z = encodeVarint(encodeZigZag(deltaZ));

                append(twkb, x, y, z);

                xShrinked = shrink(locations[i].x, precisionXY);
                yShrinked = shrink(locations[i].y, precisionXY);
                zShrinked = shrink(locations[i].z, precisionZ);

            }

            return twkb;
        }

        bytes_t makeLine(vector<PosXYZT> locations, char precisionXY, char precisionZ, char precisionT, bool bbox) {
            bytes_t twkb({0x00, 0x00});

            // Set type and precision
            setTypeAndPrecision(twkb, 0x02, precisionXY);

            // Prepare metadata header
            // bbox | size | idlist | extended dimensions | empty geom
            setMetadataBits(twkb, bbox, false, false, true, false);

            // Add extendet dimension information
            addZTDimensions(twkb, precisionZ, precisionT);

            // Add bounding box - if selected
            if (bbox) {

                int minx = shrink(locations.front().x, precisionXY);
                int maxx = minx;
                int miny = shrink(locations.front().y, precisionXY);
                int maxy = miny;
                int minz = shrink(locations.front().z, precisionZ);
                int maxz = minz;
                int mint = shrink(locations.front().t, precisionT);
                int maxt = mint;

                for (auto &location : locations) {
                    int tmpX = shrink(location.x, precisionXY);
                    int tmpY = shrink(location.y, precisionXY);
                    int tmpZ = shrink(location.z, precisionZ);
                    int tmpT = shrink(location.t, precisionT);

                    if (tmpX < minx) minx = tmpX;
                    if (tmpY < miny) miny = tmpY;
                    if (tmpZ < minz) minz = tmpZ;
                    if (tmpT < mint) mint = tmpT;
                    if (tmpX > maxx) maxx = tmpX;
                    if (tmpY > maxy) maxy = tmpY;
                    if (tmpZ > maxz) maxz = tmpZ;
                    if (tmpT > maxt) maxt = tmpT;
                }

                bytes_t bboxX = encodeVarint(encodeZigZag(minx));
                bytes_t bboxY = encodeVarint(encodeZigZag(miny));
                bytes_t bboxZ = encodeVarint(encodeZigZag(minz));
                bytes_t bboxT = encodeVarint(encodeZigZag(mint));
                bytes_t deltaX = encodeVarint(encodeZigZag(maxx - minx));
                bytes_t deltaY = encodeVarint(encodeZigZag(maxy - miny));
                bytes_t deltaZ = encodeVarint(encodeZigZag(maxz - minz));
                bytes_t deltaT = encodeVarint(encodeZigZag(maxt - mint));

                append(twkb, bboxX, deltaX, bboxY, deltaY, bboxZ, deltaZ, bboxT, deltaT);

            }


            // Set number of containing points
            bytes_t npoints = encodeVarint(locations.size());
            twkb.insert(twkb.end(), npoints.begin(), npoints.end());

            int xShrinked = shrink(locations.front().x, precisionXY);
            int yShrinked = shrink(locations.front().y, precisionXY);
            int zShrinked = shrink(locations.front().z, precisionZ);
            int tShrinked = shrink(locations.front().t, precisionT);

            bytes_t x = encodeVarint(encodeZigZag(xShrinked));
            bytes_t y = encodeVarint(encodeZigZag(yShrinked));
            bytes_t z = encodeVarint(encodeZigZag(zShrinked));
            bytes_t t = encodeVarint(encodeZigZag(tShrinked));

            append(twkb, x, y, z, t);

            for (int i = 1; i < locations.size(); i++) {

                int deltaX = shrink(locations[i].x, precisionXY) - xShrinked;
                int deltaY = shrink(locations[i].y, precisionXY) - yShrinked;
                int deltaZ = shrink(locations[i].z, precisionZ) - zShrinked;
                int deltaT = shrink(locations[i].t, precisionT) - tShrinked;

                x = encodeVarint(encodeZigZag(deltaX));
                y = encodeVarint(encodeZigZag(deltaY));
                z = encodeVarint(encodeZigZag(deltaZ));
                t = encodeVarint(encodeZigZag(deltaT));

                append(twkb, x, y, z, t);

                xShrinked = shrink(locations[i].x, precisionXY);
                yShrinked = shrink(locations[i].y, precisionXY);
                zShrinked = shrink(locations[i].z, precisionZ);
                tShrinked = shrink(locations[i].t, precisionT);

            }

            return twkb;
        }

        bytes_t makePolygon(vector<vector<PosXY>> rings, char precisionXY) {
            bytes_t twkb({0x00, 0x00});

            // Set type and precision
            setTypeAndPrecision(twkb, 0x03, precisionXY);

            // Prepare metadata header
            // bbox | size | idlist | extended dimensions | empty geom
            setMetadataBits(twkb, false, false, false, false, false);

            // Insert number of containing rings
            bytes_t nrings = encodeVarint(rings.size());
            append(twkb, nrings);

            int lastXFull = shrink(rings.front().front().x, precisionXY);
            int lastYFull = shrink(rings.front().front().y, precisionXY);

            int deltaX = lastXFull;
            int deltaY = lastYFull;

            for (size_t i = 0; i < rings.size(); i++) {
                auto &pointsInRing = rings[i];

                // Set number of containing pointsInRing
                bytes_t npoints = encodeVarint(pointsInRing.size());
                append(twkb, npoints);

                for (int j = 0; j < pointsInRing.size(); j++) {

                    int currentXFull = shrink(pointsInRing[j].x, precisionXY);
                    int currentYFull = shrink(pointsInRing[j].y, precisionXY);

                    if (!(i == 0 && j == 0)) {
                        deltaX = currentXFull - lastXFull;
                        deltaY = currentYFull - lastYFull;
                    }

                    lastXFull = currentXFull;
                    lastYFull = currentYFull;

                    bytes_t deltaXAsZigZag = encodeVarint(encodeZigZag(deltaX));
                    bytes_t deltaYAsZigZag = encodeVarint(encodeZigZag(deltaY));

                    append(twkb, deltaXAsZigZag, deltaYAsZigZag);

                }
            }


            return twkb;
        }
    };

};


#endif //TWKB_HEADER_ONLY_TWKB_H

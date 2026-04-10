from koth_ff import detect_hills, detect_features

def main():
    hills = detect_hills("/home/patrick-garrett/Data/Arabela/seminal_plasma/raw/DIA/boar/20220217_Boar-1_S3-A7_1_8655.d", mz_tolerance=8.0)
    print(hills)
    features = detect_features(
        hills)
    print(features)


if __name__ == "__main__":
    main()

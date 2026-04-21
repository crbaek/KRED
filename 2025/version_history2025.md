## 2025 update

### Summary
The 2025 update keeps the overall KRED structure unchanged at 125 monthly series, but revises several series definitions and extends the empirical-analysis sample through 2025:12. The balanced panel used in the paper is extended from 2009:06–2024:12 to 2009:06–2025:12, and `ACOGNO` is excluded from the updated balanced panel.

### Main changes

- **Balanced panel extended**
  - Empirical-analysis sample extended from **2009:06–2024:12** to **2009:06–2025:12**. :contentReference[oaicite:2]{index=2} :contentReference[oaicite:3]{index=3}

- **M1 / M2 source-side revision**
  - Raw **M1** and **M2** series are updated because the **ECOS calculation method changed**.
  - This affects the underlying money-stock data used in KRED and therefore also affects downstream transformed series based on these inputs.

- **`DTCOLNVHFNM`** is clarified from total registered motor vehicles to **private-use vehicles**.

- **`ACOGNO` definition revised**
  - `ACOGNO` is redefined in the 2025 version.
  - Previous version: proxy based on **consumer-goods imports**.
  - Updated version: proxy based on the **KOSIS Manufacturing Domestic Supply Index (consumer goods, domestic supply)**, with an explicit note that it is **not a direct new-orders series**. :contentReference[oaicite:4]{index=4} :contentReference[oaicite:5]{index=5}

- **`INVEST` construction documented more explicitly**
  - The `INVEST` series is now documented using an **explicit sum of component codes** from ECOS rather than only a high-level source label.
  - This improves transparency and reproducibility of the series construction. :contentReference[oaicite:6]{index=6}

- **`PPICMM`** is revised.
  - 2024: proxy based on **non-ferrous metal bar & basic products**
  - 2025: proxy based on the **average of basic metal products and fabricated metal products**

- **`ICSA`** description is refined.
  - 2024: `Unemployment Benefit Payment Status (Monthly)`
  - 2025: `Unemployment Benefit Claims (Monthly)`


### Additional note
- Because `ACOGNO` is now treated differently and has incomplete coverage under the updated setup, it is excluded from the 2025 balanced panel used for factor analysis. 



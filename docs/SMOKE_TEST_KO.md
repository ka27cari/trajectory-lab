# Trajectory Lab v1.2.0: 빠른 화면 확인

R 또는 RStudio에서 프로젝트 폴더를 연 뒤 `source("install.R")`, `shiny::runApp(".")`를 실행한다.

1. **Correlations:** 기본 시뮬레이션에서 상관행렬이 보이고 `MEAN |r|`가 `N/A`가 아닌지 확인한다. `Spearman`, 시간 가중치를 선택해도 실행되는지 확인한다.
2. **PCA:** `Feature-by-feature PCA`에서 `Previous feature`, `Add next feature`를 누르면 PC1/PC2 점 경로와 오른쪽 누적 loading 경로가 바뀌는지 확인한다. 두 계산 방식(`Refit`, `Fixed`)을 바꿔본다.
3. **Mean & variance:** 기본값(offset=0, spread=1)에서는 실제 평균 이동량이 0이고 분산비가 1인지 확인한다. `Requested mean offset=0.1`, `Spread multiplier=1`에서는 평균만 바뀐다(값이 0과 1에 닿으면 clipping 예외).
4. **분산 압축:** offset=0, spread=0.5를 설정하면 clipping이 없는 feature의 분산비가 0.25에 가까운지 확인한다. 분산비는 SD 비율과 다르다.
5. **적용 스위치:** `Apply this experiment...`를 켜고 PCA/Correlations 탭을 열면 분석 결과가 바뀌는지, 끄면 원래 결과로 돌아오는지 확인한다.
6. **글자 겹침:** 브라우저 너비를 줄이고 time colorbar 눈금/heatmap 축 이름이 겹치지 않는지, feature를 30개까지 늘려 다시 확인한다.
7. **Distance maps (Euclidean):** 기본값에서 `Recalculate four maps`를 눌러 2×2 비교 화면, 품질표, 세부 그래프, 거리 히트맵이 보이는지 확인한다. 시간점 표본 수(기본 120)와 경고문을 확인한다.
8. **다섯 가지 거리:** Manhattan/Canberra/상관거리로 바꾸고 재계산한다. Bray-Curtis는 `Standardize`를 끄고 재계산한다. 나머지 조건이 같으면 PCA 좌표는 거리 종류만 바꿔도 변하지 않는다.
9. **매핑 대상:** `Objects being mapped = Features`를 선택해 각 CpG 프로필이 점 하나가 되는지 확인한다. 표본 수를 줄이거나 실험 스위치를 켠 뒤에는 반드시 `Recalculate`를 누른다.
10. **거리 실패 처리:** 시계열 행을 동일하게 만든 CSV 등으로 0 거리의 객체가 있을 때 NMDS/Sammon만 `Unavailable`로 표시되고 PCA/PCoA는 계속 그려지는지 확인한다.
11. **글자 겹침:** PC1/PC2 그래프, 2×2 비교, 거리 히트맵, 상관 히트맵을 모바일 화면과 30개 feature로 확인한다.
12. **Export:** `PCA step · CSV`와 `PCA step · PNG`, `Mean & variance · CSV`, `Active window · CSV` 및 거리 좌표/품질표/행렬/PNG가 다운로드되는지 확인한다.

기본 테스트 `Rscript tests/run_tests.R`가 성공한 뒤에 배포한다. 이 테스트는 계산 및 파싱 검증용이고, 웹 브라우저 렌더링을 자동 검증하지는 않는다.

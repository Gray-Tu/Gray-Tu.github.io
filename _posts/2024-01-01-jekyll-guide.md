---
title: "部落格使用指南"
categories: [筆記]
tags: [Jekyll, Markdown, GitHub Pages]
---

這篇文章說明如何使用這個 Jekyll 部落格系統。

## 新增文章

在 `_posts/` 資料夾中新建一個 Markdown 檔案：

```
_posts/2024-01-15-my-new-post.md
```

## Front Matter

每篇文章開頭需要加上 YAML front matter：

```yaml
---
title: "文章標題"
categories: [分類]
tags: [標籤1, 標籤2]
---
```

## Markdown 語法

支援完整的 Markdown 語法，包含：

- **粗體** 和 *斜體*
- [連結](https://example.com)
- 程式碼區塊
- 表格、引用、清單等

### 程式碼範例

```python
import pandas as pd

df = pd.read_csv("data.csv")
print(df.head())
```

### 引用

> 這是一段引用文字，適合用來標記重要的筆記或參考資料。

## 圖片

將圖片放在 `assets/images/` 資料夾，然後在文章中引用：

```markdown
![圖片說明](/assets/images/your-image.png)
```

## 分類與標籤

- **分類 (categories)**：用於文章的主要分類，例如「生物資訊」、「筆記」
- **標籤 (tags)**：更細緻的標記，例如「Python」、「R」、「10x」

## 發布

將變更推上 GitHub 即可：

```bash
git add _posts/your-new-post.md
git commit -m "新增文章：文章標題"
git push
```

GitHub Pages 會自動建置並發布你的網站。

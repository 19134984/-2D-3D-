-- Gives every example in the PDF a fixed LaTeX label.
--
-- link_examples.sh wraps each example it includes into the book in
--
--   ::: {.olb-example label="ex:<path>"}
--   {{< include <path>.qmd >}}
--   :::
--
-- This filter puts \label{ex:<path>} right after the first heading of the
-- included page, so that \ref and \hyperref -- used by the overview table
-- that example_table.py writes -- resolve to the section of the example, and
-- drops the div. The heading keeps the id Quarto gave it.

function Div(div)
  if not div.classes:includes("olb-example") then
    return nil
  end
  local label = div.attributes["label"]
  local blocks = div.content
  for i, block in ipairs(blocks) do
    if block.t == "Header" then
      blocks:insert(i + 1, pandoc.RawBlock("latex", "\\label{" .. label .. "}"))
      return blocks
    end
  end
  quarto.log.warning("example_labels: no heading in the page of " .. label)
  return blocks
end

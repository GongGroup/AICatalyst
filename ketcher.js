function getInput() {
  var url = location.search;
  var theInput = new Object();
  if (url.indexOf("?") != -1) {
    var str = url.slice(1);
    strs = str.split("&");
    for (var i = 0; i < strs.length; i++) {
      theInput[strs[i].split("=")[0]] = decodeURI(strs[i].split("=")[1]);
    }
  }
  return theInput;
}

function handleInput() {
  var theInput = getInput();
  if (Object.keys(theInput).length > 0) {
    name2structure(theInput["input"]);
  }
}

function getKetcher() {
  var ketcherFrame = document.getElementById("ifKetcher");
  var ketcher;
  if ("contentDocument" in ketcherFrame) {
    ketcher = ketcherFrame.contentWindow.ketcher;
  } // IE7
  else {
    ketcher = document.frames["ifKetcher"].window.ketcher;
  }

  return ketcher;
}

function name2structure(chemicalName) {
  $.ajax({
    url: "../chemical/record.json",
    type: "GET",
    dataType: "json",
    success: function (data) {
      var emol = data[chemicalName].replace("Ok.", "").replace(/\r\r/g, "");
      // console.log(emol);
      var ketcher = getKetcher();
      ketcher.setMolecule(emol);
    },
  });
}

function structure2smile() {
  var ketcher = getKetcher();
  var iename = document.getElementById("iename");
  p = ketcher.getSmiles();
  p.then((name) => (iename.innerHTML = name));
}
